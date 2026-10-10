//! Hand-rolled tokenizer + recursive-descent parser for the bounded
//! Abaqus scripting subset accepted by `abaqus2solvix` (issue #35).
//!
//! Documented boundary — this is NOT general Python:
//! * expression statements and simple `name = expr` assignments;
//! * attribute chains (`mdb.models['M'].parts['P']...`), string/number
//!   subscriptions, method calls with keyword arguments (positional
//!   arguments parse; the translator warns when it must ignore them);
//! * numbers, strings, lists, tuples, unary minus, `+ - * /`;
//! * `#` comments, bracketed multi-line continuations, backslash
//!   continuations, `;` statement separators, `import` lines;
//! * no control flow, classes, functions or slicing: `if`/`for`/
//!   `while`/`def`/... blocks are skipped whole (indentation-tracked)
//!   and reported by the translator as unsupported, so nothing inside
//!   a loop is silently mistranslated.

use std::collections::HashMap;

/// Parse error: 1-based source line + message.
#[derive(Debug, Clone, PartialEq)]
pub struct PyError {
    pub line: usize,
    pub message: String,
}

fn perr<T>(line: usize, message: impl Into<String>) -> Result<T, PyError> {
    Err(PyError {
        line,
        message: message.into(),
    })
}

#[derive(Debug, Clone)]
pub struct Script {
    pub stmts: Vec<Stmt>,
}

#[derive(Debug, Clone)]
pub struct Stmt {
    pub line: usize,
    pub kind: StmtKind,
}

#[derive(Debug, Clone)]
pub enum StmtKind {
    /// `import x` / `from abaqus import *` — recognized, ignored.
    Import,
    /// `pass` — no-op.
    Noop,
    /// A control-flow block header: the whole indented block was
    /// skipped. The translator reports it as unsupported.
    Block { keyword: String },
    /// `target = value` (simple name binding only).
    Assign { target: String, value: Expr },
    /// Bare expression statement (the Abaqus API calls).
    Expr(Expr),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Num(f64),
    Str(String),
    Name(String),
    List(Vec<Expr>),
    Tuple(Vec<Expr>),
    Attr(Box<Expr>, String),
    Index(Box<Expr>, Box<Expr>),
    Call {
        func: Box<Expr>,
        pos: Vec<Expr>,
        kwargs: Vec<(String, Expr)>,
    },
    Bin(Box<Expr>, char, Box<Expr>),
    Neg(Box<Expr>),
}

/// Positional + keyword call arguments.
pub type CallArgs = (Vec<Expr>, Vec<(String, Expr)>);

/// Plain-data Python value (everything the evaluator can compute).
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    None,
    Bool(bool),
    Num(f64),
    Str(String),
    List(Vec<Value>),
    Tuple(Vec<Value>),
}

impl Value {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Num(x) => Some(*x),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }
}

/// Evaluate a pure-data expression (numbers, strings, lists, tuples,
/// arithmetic, bound variables). API object expressions
/// (attribute/index/call chains) are rejected.
pub fn eval_expr(e: &Expr, env: &HashMap<String, Value>) -> Result<Value, String> {
    match e {
        Expr::Num(x) => Ok(Value::Num(*x)),
        Expr::Str(s) => Ok(Value::Str(s.clone())),
        Expr::Name(n) => match n.as_str() {
            "True" => Ok(Value::Bool(true)),
            "False" => Ok(Value::Bool(false)),
            "None" => Ok(Value::None),
            _ => env
                .get(n)
                .cloned()
                .ok_or_else(|| format!("unbound name '{n}'")),
        },
        Expr::List(items) => {
            let mut out = Vec::with_capacity(items.len());
            for it in items {
                out.push(eval_expr(it, env)?);
            }
            Ok(Value::List(out))
        }
        Expr::Tuple(items) => {
            let mut out = Vec::with_capacity(items.len());
            for it in items {
                out.push(eval_expr(it, env)?);
            }
            Ok(Value::Tuple(out))
        }
        Expr::Bin(a, op, b) => {
            let x = eval_expr(a, env)?.as_f64().ok_or("non-numeric operand")?;
            let y = eval_expr(b, env)?.as_f64().ok_or("non-numeric operand")?;
            match op {
                '+' => Ok(Value::Num(x + y)),
                '-' => Ok(Value::Num(x - y)),
                '*' => Ok(Value::Num(x * y)),
                '/' if y == 0.0 => Err("division by zero".into()),
                '/' => Ok(Value::Num(x / y)),
                c => Err(format!("unsupported operator '{c}'")),
            }
        }
        Expr::Neg(a) => Ok(Value::Num(
            -eval_expr(a, env)?.as_f64().ok_or("non-numeric operand")?,
        )),
        Expr::Attr { .. } | Expr::Index { .. } | Expr::Call { .. } => {
            Err("API expression used where a data value is expected".into())
        }
    }
}

// ---------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Ident(String),
    Num(f64),
    Str(String),
    Punct(char),
    Newline,
}

#[derive(Debug, Clone)]
struct Token {
    tok: Tok,
    line: usize,
    col: usize,
}

const PUNCTS: &[char] = &[
    '(', ')', '[', ']', '{', '}', ',', '.', '=', ':', ';', '+', '-', '*', '/',
];

#[allow(unused_assignments)] // col=1 after the final newline is a harmless dead store
fn tokenize(src: &str) -> Result<Vec<Token>, PyError> {
    let ch: Vec<char> = src.chars().collect();
    let mut toks = Vec::new();
    let (mut line, mut col) = (1usize, 1usize);
    let mut depth = 0i32;
    let mut i = 0usize;
    let n = ch.len();
    while i < n {
        let c = ch[i];
        match c {
            '\n' => {
                i += 1;
                if depth == 0 {
                    toks.push(Token {
                        tok: Tok::Newline,
                        line,
                        col,
                    });
                }
                line += 1;
                col = 1;
            }
            ' ' | '\t' | '\r' => {
                i += 1;
                col += 1;
            }
            '#' => {
                while i < n && ch[i] != '\n' {
                    i += 1;
                    col += 1;
                }
            }
            '\\' => {
                i += 1;
                col += 1;
                if i < n && ch[i] == '\r' {
                    i += 1;
                    col += 1;
                }
                if i < n && ch[i] == '\n' {
                    i += 1;
                    line += 1;
                    col = 1;
                } else {
                    return perr(line, "backslash continuation must end a line");
                }
            }
            '\'' | '"' => {
                let quote = c;
                let (mut j, mut s) = (i + 1, String::new());
                while j < n && ch[j] != quote {
                    if ch[j] == '\n' {
                        return perr(line, "unterminated string");
                    }
                    if ch[j] == '\\' {
                        j += 1;
                        if j >= n {
                            return perr(line, "unterminated string escape");
                        }
                        s.push(match ch[j] {
                            'n' => '\n',
                            't' => '\t',
                            'r' => '\r',
                            '\\' => '\\',
                            '\'' => '\'',
                            '"' => '"',
                            other => other,
                        });
                    } else {
                        s.push(ch[j]);
                    }
                    j += 1;
                }
                if j >= n {
                    return perr(line, "unterminated string");
                }
                col += j + 1 - i;
                i = j + 1;
                toks.push(Token {
                    tok: Tok::Str(s),
                    line,
                    col,
                });
            }
            c if c.is_ascii_digit() || (c == '.' && i + 1 < n && ch[i + 1].is_ascii_digit()) => {
                let start = i;
                let mut seen_dot = false;
                while i < n && (ch[i].is_ascii_digit() || (ch[i] == '.' && !seen_dot)) {
                    if ch[i] == '.' {
                        seen_dot = true;
                    }
                    i += 1;
                }
                if i < n && (ch[i] == 'e' || ch[i] == 'E') {
                    i += 1;
                    if i < n && (ch[i] == '+' || ch[i] == '-') {
                        i += 1;
                    }
                    let d0 = i;
                    while i < n && ch[i].is_ascii_digit() {
                        i += 1;
                    }
                    if i == d0 {
                        return perr(line, "malformed number exponent");
                    }
                }
                let text: String = ch[start..i].iter().collect();
                let x = text.parse::<f64>().map_err(|_| PyError {
                    line,
                    message: format!("malformed number '{text}'"),
                })?;
                col += i - start;
                toks.push(Token {
                    tok: Tok::Num(x),
                    line,
                    col,
                });
            }
            c if c.is_ascii_alphabetic() || c == '_' => {
                let start = i;
                while i < n && (ch[i].is_ascii_alphanumeric() || ch[i] == '_') {
                    i += 1;
                }
                let name: String = ch[start..i].iter().collect();
                col += i - start;
                toks.push(Token {
                    tok: Tok::Ident(name),
                    line,
                    col,
                });
            }
            '(' | '[' | '{' => {
                depth += 1;
                toks.push(Token {
                    tok: Tok::Punct(c),
                    line,
                    col,
                });
                i += 1;
                col += 1;
            }
            ')' | ']' | '}' => {
                depth -= 1;
                if depth < 0 {
                    return perr(line, "unbalanced closing bracket");
                }
                toks.push(Token {
                    tok: Tok::Punct(c),
                    line,
                    col,
                });
                i += 1;
                col += 1;
            }
            c if PUNCTS.contains(&c) => {
                toks.push(Token {
                    tok: Tok::Punct(c),
                    line,
                    col,
                });
                i += 1;
                col += 1;
            }
            _ => return perr(line, format!("unexpected character '{c}'")),
        }
    }
    toks.push(Token {
        tok: Tok::Newline,
        line,
        col,
    });
    Ok(toks)
}

/// Split the token stream into top-level statement groups.
fn statements(toks: Vec<Token>) -> Vec<Vec<Token>> {
    let mut groups = Vec::new();
    let mut cur: Vec<Token> = Vec::new();
    for t in toks {
        match t.tok {
            Tok::Newline | Tok::Punct(';') => {
                if !cur.is_empty() {
                    groups.push(std::mem::take(&mut cur));
                }
            }
            _ => cur.push(t),
        }
    }
    if !cur.is_empty() {
        groups.push(cur);
    }
    groups
}

// ---------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------

const BLOCK_KEYWORDS: &[&str] = &[
    "if", "elif", "else", "for", "while", "def", "class", "with", "try", "except", "finally",
];

pub fn parse_script(src: &str) -> Result<Script, PyError> {
    let groups = statements(tokenize(src)?);
    let mut stmts = Vec::new();
    let mut i = 0;
    while i < groups.len() {
        let first = &groups[i][0];
        if let Tok::Ident(kw) = &first.tok {
            if BLOCK_KEYWORDS.contains(&kw.as_str()) {
                stmts.push(Stmt {
                    line: first.line,
                    kind: StmtKind::Block {
                        keyword: kw.clone(),
                    },
                });
                let ind = first.col;
                i += 1;
                while i < groups.len() && groups[i][0].col > ind {
                    i += 1;
                }
                continue;
            }
        }
        stmts.push(parse_stmt(&groups[i])?);
        i += 1;
    }
    Ok(Script { stmts })
}

fn parse_stmt(toks: &[Token]) -> Result<Stmt, PyError> {
    let line = toks[0].line;
    if let Tok::Ident(w) = &toks[0].tok {
        match w.as_str() {
            "import" | "from" => {
                return Ok(Stmt {
                    line,
                    kind: StmtKind::Import,
                })
            }
            "pass" => {
                return Ok(Stmt {
                    line,
                    kind: StmtKind::Noop,
                })
            }
            "del" | "global" | "nonlocal" | "assert" | "raise" | "return" | "yield" | "lambda" => {
                return perr(line, format!("'{w}' is outside the translatable subset"));
            }
            _ => {}
        }
    }
    let mut p = Parser { toks, pos: 0 };
    let e = p.parse_expr()?;
    if p.peek_punct('=') {
        let Expr::Name(target) = &e else {
            return perr(line, "only simple `name = value` assignments are supported");
        };
        p.pos += 1;
        let v = p.parse_expr()?;
        p.expect_end(line)?;
        Ok(Stmt {
            line,
            kind: StmtKind::Assign {
                target: target.clone(),
                value: v,
            },
        })
    } else {
        p.expect_end(line)?;
        Ok(Stmt {
            line,
            kind: StmtKind::Expr(e),
        })
    }
}

struct Parser<'a> {
    toks: &'a [Token],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|t| &t.tok)
    }

    fn peek_at(&self, k: usize) -> Option<&Tok> {
        self.toks.get(self.pos + k).map(|t| &t.tok)
    }

    fn line(&self) -> usize {
        self.toks
            .get(self.pos)
            .or_else(|| self.toks.last())
            .map(|t| t.line)
            .unwrap_or(0)
    }

    fn bump_tok(&mut self) -> Option<Tok> {
        let t = self.toks.get(self.pos).map(|x| x.tok.clone());
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn peek_punct(&self, c: char) -> bool {
        matches!(self.peek(), Some(Tok::Punct(p)) if *p == c)
    }

    fn eat_punct(&mut self, c: char) -> bool {
        if self.peek_punct(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect_punct(&mut self, c: char) -> Result<(), PyError> {
        if self.eat_punct(c) {
            Ok(())
        } else {
            perr(self.line(), format!("expected '{c}'"))
        }
    }

    fn expect_end(&self, line: usize) -> Result<(), PyError> {
        if self.pos >= self.toks.len() {
            Ok(())
        } else {
            perr(line, "unexpected trailing tokens in statement")
        }
    }

    fn parse_expr(&mut self) -> Result<Expr, PyError> {
        let mut a = self.parse_term()?;
        while let Some(Tok::Punct(p)) = self.peek().cloned() {
            if p != '+' && p != '-' {
                break;
            }
            self.pos += 1;
            let b = self.parse_term()?;
            a = Expr::Bin(Box::new(a), p, Box::new(b));
        }
        Ok(a)
    }

    fn parse_term(&mut self) -> Result<Expr, PyError> {
        let mut a = self.parse_unary()?;
        while let Some(Tok::Punct(p)) = self.peek().cloned() {
            if p != '*' && p != '/' {
                break;
            }
            self.pos += 1;
            let b = self.parse_unary()?;
            a = Expr::Bin(Box::new(a), p, Box::new(b));
        }
        Ok(a)
    }

    fn parse_unary(&mut self) -> Result<Expr, PyError> {
        if self.eat_punct('-') {
            Ok(Expr::Neg(Box::new(self.parse_unary()?)))
        } else {
            self.parse_postfix()
        }
    }

    fn parse_postfix(&mut self) -> Result<Expr, PyError> {
        let mut e = self.parse_atom()?;
        loop {
            if self.eat_punct('.') {
                let name = match self.bump_tok() {
                    Some(Tok::Ident(n)) => n,
                    _ => return perr(self.line(), "expected attribute name after '.'"),
                };
                e = Expr::Attr(Box::new(e), name);
            } else if self.eat_punct('[') {
                let idx = self.parse_expr()?;
                self.expect_punct(']')?;
                e = Expr::Index(Box::new(e), Box::new(idx));
            } else if self.peek_punct('(') {
                let (pos, kwargs) = self.parse_args()?;
                e = Expr::Call {
                    func: Box::new(e),
                    pos,
                    kwargs,
                };
            } else {
                break;
            }
        }
        Ok(e)
    }

    fn parse_atom(&mut self) -> Result<Expr, PyError> {
        match self.bump_tok() {
            Some(Tok::Num(x)) => Ok(Expr::Num(x)),
            Some(Tok::Str(s)) => Ok(Expr::Str(s)),
            Some(Tok::Ident(n)) => Ok(Expr::Name(n)),
            Some(Tok::Punct('(')) => self.parse_paren(),
            Some(Tok::Punct('[')) => self.parse_list(),
            _ => perr(self.line(), "expected a value"),
        }
    }

    fn parse_paren(&mut self) -> Result<Expr, PyError> {
        if self.eat_punct(')') {
            return Ok(Expr::Tuple(vec![]));
        }
        let mut items = vec![self.parse_expr()?];
        let mut trailing = false;
        while self.eat_punct(',') {
            trailing = true;
            if self.peek_punct(')') {
                break;
            }
            items.push(self.parse_expr()?);
            trailing = false;
        }
        self.expect_punct(')')?;
        if items.len() == 1 && !trailing {
            Ok(items.pop().unwrap())
        } else {
            Ok(Expr::Tuple(items))
        }
    }

    fn parse_list(&mut self) -> Result<Expr, PyError> {
        if self.eat_punct(']') {
            return Ok(Expr::List(vec![]));
        }
        let mut items = vec![self.parse_expr()?];
        while self.eat_punct(',') {
            if self.peek_punct(']') {
                break;
            }
            items.push(self.parse_expr()?);
        }
        self.expect_punct(']')?;
        Ok(Expr::List(items))
    }

    fn parse_args(&mut self) -> Result<CallArgs, PyError> {
        self.expect_punct('(')?;
        let mut pos = Vec::new();
        let mut kwargs = Vec::new();
        if self.eat_punct(')') {
            return Ok((pos, kwargs));
        }
        loop {
            if matches!(self.peek(), Some(Tok::Ident(_)))
                && matches!(self.peek_at(1), Some(Tok::Punct('=')))
            {
                let Some(Tok::Ident(name)) = self.bump_tok() else {
                    break;
                };
                self.pos += 1; // '='
                let v = self.parse_expr()?;
                kwargs.push((name, v));
            } else if self.peek_punct('*') {
                return perr(self.line(), "*args/**kwargs are not supported");
            } else {
                pos.push(self.parse_expr()?);
            }
            if !self.eat_punct(',') {
                break;
            }
            if self.peek_punct(')') {
                break;
            }
        }
        self.expect_punct(')')?;
        Ok((pos, kwargs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_calls_kwargs_and_tuples() {
        let s = parse_script(
            "from abaqus import *\n\
             s.rectangle(point1=(0.0, -0.05), point2=(1.0, 0.05))\n\
             mdb.models['M'].materials['Steel'].Elastic(table=((210e9, 0.3), ))\n",
        )
        .unwrap();
        assert_eq!(s.stmts.len(), 3);
        assert!(matches!(s.stmts[0].kind, StmtKind::Import));
        assert!(matches!(s.stmts[1].kind, StmtKind::Expr(_)));
    }

    #[test]
    fn bracket_continuation_joins_lines() {
        let s = parse_script(
            "p = mdb.models['Model-1'].Part(name='Beam',\n\
             dimensionality=TWO_D_PLANAR,\n\
             type=DEFORMABLE_BODY)\n\
             x = 5\n",
        )
        .unwrap();
        assert_eq!(s.stmts.len(), 2);
        let StmtKind::Assign { value, .. } = &s.stmts[0].kind else {
            panic!("expected assignment");
        };
        let expected = Expr::Call {
            func: Box::new(Expr::Attr(
                Box::new(Expr::Index(
                    Box::new(Expr::Attr(
                        Box::new(Expr::Name("mdb".into())),
                        "models".into(),
                    )),
                    Box::new(Expr::Str("Model-1".into())),
                )),
                "Part".into(),
            )),
            pos: vec![],
            kwargs: vec![
                ("name".into(), Expr::Str("Beam".into())),
                ("dimensionality".into(), Expr::Name("TWO_D_PLANAR".into())),
                ("type".into(), Expr::Name("DEFORMABLE_BODY".into())),
            ],
        };
        // Index key is the *string* 'Model-1'; the func chain mirrors it.
        assert_eq!(*value, expected);
    }

    #[test]
    fn control_flow_blocks_are_skipped_whole() {
        let s = parse_script(
            "for i in range(3):\n    s.Line(point1=(i, 0), point2=(i, 1))\n    x = i\n\
             y = 2\n",
        )
        .unwrap();
        assert_eq!(s.stmts.len(), 2);
        assert!(matches!(
            &s.stmts[0].kind,
            StmtKind::Block { keyword } if keyword == "for"
        ));
        assert!(matches!(&s.stmts[1].kind, StmtKind::Assign { target, .. } if target == "y"));
    }

    #[test]
    fn evaluates_arithmetic_and_env() {
        let mut env = HashMap::new();
        env.insert("L".to_string(), Value::Num(2.0));
        let e = parse_script("x = L / 4 + 1\n").unwrap();
        let StmtKind::Assign { value, .. } = e.stmts[0].kind.clone() else {
            panic!();
        };
        assert_eq!(eval_expr(&value, &env).unwrap(), Value::Num(1.5));
    }

    #[test]
    fn tuple_trailing_comma_semantics() {
        // ((7850.0), ) is a 1-tuple whose element is a bare number.
        let s = parse_script("t = ((7850.0), )\n").unwrap();
        let StmtKind::Assign { value, .. } = s.stmts[0].kind.clone() else {
            panic!()
        };
        let v = eval_expr(&value, &HashMap::new()).unwrap();
        assert_eq!(v, Value::Tuple(vec![Value::Num(7850.0)]));
    }
}
