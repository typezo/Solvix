//! Binary `.solvix` output database (issue #32).
//!
//! Container layout — all values little-endian, no external deps:
//!
//! ```text
//! magic "SLVXODB1"          8 bytes
//! u32 version (= 1)
//! model: u32 n_nodes, then per node (u32 id, f64 x, f64 y, f64 z);
//!        u32 n_elements, then per element (u32 id, u8 kind,
//!        u32 conn_len, conn_len x u32 node ids)
//! frames (appended by OdbWriter::append_frame):
//!   per frame: u32 flags (bit0 temperatures, bit1 velocities),
//!   f64 time, f64 load_fraction, then sized arrays:
//!   displacements (3 x f64 per node), element von Mises (f64 each),
//!   nodal von Mises (f64 each), temperatures (f64 each, bit0),
//!   velocities (3 x f64 per node, bit1)
//! history region (at finish): |u| per node per frame, node-major
//! index (at finish): u32 n_frames, u64 frame offset each,
//!   u32 flags union, u32 n_nodes, u64 history offset per node
//! trailer: u64 index offset — the LAST 8 bytes of the file
//! ```
//!
//! Random access reads seek the trailer, then the index, then frames.

use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use solvix_core::{ElemKind, Mesh};

const MAGIC: &[u8; 8] = b"SLVXODB1";
const VERSION: u32 = 1;
const FLAG_TEMPERATURES: u32 = 1;
const FLAG_VELOCITIES: u32 = 2;

fn kind_code(k: ElemKind) -> u8 {
    match k {
        ElemKind::Cps3 => 1,
        ElemKind::Cps4 => 2,
        ElemKind::Cps8 => 3,
        ElemKind::Cpe4 => 4,
        ElemKind::Cpe8 => 5,
        ElemKind::C3d8 => 6,
        ElemKind::C3d4 => 7,
        ElemKind::B31 => 8,
    }
}

fn kind_from_code(code: u8) -> io::Result<ElemKind> {
    match code {
        1 => Ok(ElemKind::Cps3),
        2 => Ok(ElemKind::Cps4),
        3 => Ok(ElemKind::Cps8),
        4 => Ok(ElemKind::Cpe4),
        5 => Ok(ElemKind::Cpe8),
        6 => Ok(ElemKind::C3d8),
        7 => Ok(ElemKind::C3d4),
        8 => Ok(ElemKind::B31),
        other => Err(invalid(format!("unknown element kind code {other}"))),
    }
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// One output frame. Array lengths follow the model: displacements and
/// nodal fields parallel `mesh.nodes`, von Mises parallels
/// `mesh.elements`; the optionals are present for thermal (temperatures)
/// and dynamic (velocities) jobs.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OdbFrame {
    pub time: f64,
    /// 0..1 for NL increments and normalized dynamic time, 1.0 for a
    /// static frame.
    pub load_fraction: f64,
    pub displacements: Vec<[f64; 3]>,
    pub von_mises: Vec<f64>,
    pub nodal_von_mises: Vec<f64>,
    pub temperatures: Option<Vec<f64>>,
    pub velocities: Option<Vec<[f64; 3]>>,
}

/// Incremental `.solvix` writer: writes the model once at [`create`],
/// then appends frames; [`finish`] writes the history region, the
/// offset index and the trailer. Files are complete only after
/// `finish` (readers need the trailer).
// ponytail: history buffered in memory (8 B * nodes * frames) — stream
// it per-node when transient jobs get big enough to matter.
pub struct OdbWriter {
    file: BufWriter<std::fs::File>,
    pos: u64,
    frame_offsets: Vec<u64>,
    /// Per-frame, per-node |u| — transposed to node-major in `finish`.
    history: Vec<Vec<f64>>,
    flags_union: u32,
}

impl OdbWriter {
    /// Create (truncate) `path` and write the mesh as the model index.
    pub fn create(path: &Path, mesh: &Mesh) -> io::Result<Self> {
        let mut buf = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&VERSION.to_le_bytes());
        put_u32(&mut buf, mesh.nodes.len() as u32);
        for n in &mesh.nodes {
            put_u32(&mut buf, n.id);
            for &c in &n.coords {
                put_f64(&mut buf, c);
            }
        }
        put_u32(&mut buf, mesh.elements.len() as u32);
        for e in &mesh.elements {
            put_u32(&mut buf, e.id);
            buf.push(kind_code(e.kind));
            put_u32(&mut buf, e.connectivity.len() as u32);
            for &nid in &e.connectivity {
                put_u32(&mut buf, nid);
            }
        }
        let mut file = BufWriter::new(std::fs::File::create(path)?);
        file.write_all(&buf)?;
        Ok(Self {
            pos: buf.len() as u64,
            file,
            frame_offsets: Vec::new(),
            history: Vec::new(),
            flags_union: 0,
        })
    }

    /// Append one frame, recording its offset and per-node |u| history.
    pub fn append_frame(&mut self, frame: &OdbFrame) -> io::Result<()> {
        let mut flags = 0u32;
        if frame.temperatures.is_some() {
            flags |= FLAG_TEMPERATURES;
        }
        if frame.velocities.is_some() {
            flags |= FLAG_VELOCITIES;
        }
        let mut buf = Vec::new();
        put_u32(&mut buf, flags);
        put_f64(&mut buf, frame.time);
        put_f64(&mut buf, frame.load_fraction);
        put_u32(&mut buf, frame.displacements.len() as u32);
        for d in &frame.displacements {
            for &c in d {
                put_f64(&mut buf, c);
            }
        }
        put_u32(&mut buf, frame.von_mises.len() as u32);
        for &v in &frame.von_mises {
            put_f64(&mut buf, v);
        }
        put_u32(&mut buf, frame.nodal_von_mises.len() as u32);
        for &v in &frame.nodal_von_mises {
            put_f64(&mut buf, v);
        }
        if let Some(temps) = &frame.temperatures {
            put_u32(&mut buf, temps.len() as u32);
            for &t in temps {
                put_f64(&mut buf, t);
            }
        }
        if let Some(vels) = &frame.velocities {
            put_u32(&mut buf, vels.len() as u32);
            for v in vels {
                for &c in v {
                    put_f64(&mut buf, c);
                }
            }
        }
        self.frame_offsets.push(self.pos);
        self.history.push(
            frame
                .displacements
                .iter()
                .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
                .collect(),
        );
        self.flags_union |= flags;
        self.file.write_all(&buf)?;
        self.pos += buf.len() as u64;
        Ok(())
    }

    /// Write the history region, offset index and trailer; consume the
    /// writer. The file is readable by [`Odb`] afterwards.
    pub fn finish(mut self) -> io::Result<()> {
        let n_frames = self.frame_offsets.len() as u32;
        let n_nodes = self.history.first().map(|h| h.len()).unwrap_or(0) as u32;
        // History region, node-major: node i's full |u| series.
        let mut hist = Vec::new();
        let mut hist_offsets = Vec::with_capacity(n_nodes as usize);
        for i in 0..n_nodes as usize {
            hist_offsets.push(self.pos + (hist.len() as u64));
            for per_frame in &self.history {
                put_f64(&mut hist, per_frame.get(i).copied().unwrap_or(0.0));
            }
        }
        self.file.write_all(&hist)?;
        self.pos += hist.len() as u64;
        let index_offset = self.pos;
        let mut idx = Vec::new();
        put_u32(&mut idx, n_frames);
        for &off in &self.frame_offsets {
            put_u64(&mut idx, off);
        }
        put_u32(&mut idx, self.flags_union);
        put_u32(&mut idx, n_nodes);
        for &off in &hist_offsets {
            put_u64(&mut idx, off);
        }
        self.file.write_all(&idx)?;
        // Trailer: the index offset is the LAST 8 bytes of the file.
        let mut trailer = Vec::with_capacity(8);
        put_u64(&mut trailer, index_offset);
        self.file.write_all(&trailer)?;
        self.file.flush()?;
        Ok(())
    }
}

fn put_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_f64(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn read_exact_bytes(f: &mut impl Read, n: usize) -> io::Result<Vec<u8>> {
    let mut v = vec![0u8; n];
    f.read_exact(&mut v)?;
    Ok(v)
}

fn read_u32(f: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    f.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64(f: &mut impl Read) -> io::Result<u64> {
    let mut b = [0u8; 8];
    f.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_u8(f: &mut impl Read) -> io::Result<u8> {
    let mut b = [0u8; 1];
    f.read_exact(&mut b)?;
    Ok(b[0])
}

fn read_f64(f: &mut impl Read) -> io::Result<f64> {
    let mut b = [0u8; 8];
    f.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

fn read_f64s(f: &mut impl Read, n: usize) -> io::Result<Vec<f64>> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(read_f64(f)?);
    }
    Ok(v)
}

fn read_triples(f: &mut impl Read, n: usize) -> io::Result<Vec<[f64; 3]>> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push([read_f64(f)?, read_f64(f)?, read_f64(f)?]);
    }
    Ok(v)
}

/// Random-access reader for a finished `.solvix` database: seeks the
/// end-of-file trailer, loads the offset index and the model once, then
/// serves frames and history by seek.
#[derive(Debug)]
pub struct Odb {
    file: std::fs::File,
    model: Mesh,
    frame_offsets: Vec<u64>,
    history_offsets: Vec<u64>,
    flags: u32,
    bytes: u64,
}

impl Odb {
    /// Open and index a database file. `InvalidData` for a wrong magic,
    /// unsupported version or malformed index; other `io::Error`s pass
    /// through.
    pub fn open(path: &Path) -> io::Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let bytes = file.metadata()?.len();
        if bytes < MAGIC.len() as u64 + 4 + 8 {
            return Err(invalid("file too small to be a .solvix database"));
        }
        if read_exact_bytes(&mut file, 8)?.as_slice() != MAGIC {
            return Err(invalid("bad magic: not a .solvix database"));
        }
        let version = read_u32(&mut file)?;
        if version != VERSION {
            return Err(invalid(format!("unsupported .solvix version {version}")));
        }
        let model = read_model(&mut file)?;
        // Trailer: last 8 bytes point at the index.
        file.seek(SeekFrom::End(-8))?;
        let index_offset = read_u64(&mut file)?;
        if index_offset >= bytes {
            return Err(invalid("index offset beyond end of file"));
        }
        file.seek(SeekFrom::Start(index_offset))?;
        let n_frames = read_u32(&mut file)? as usize;
        let mut frame_offsets = Vec::with_capacity(n_frames);
        for _ in 0..n_frames {
            let off = read_u64(&mut file)?;
            if off >= bytes {
                return Err(invalid("frame offset beyond end of file"));
            }
            frame_offsets.push(off);
        }
        let flags = read_u32(&mut file)?;
        let n_nodes = read_u32(&mut file)? as usize;
        if n_nodes != model.nodes.len() {
            return Err(invalid("history node count does not match model"));
        }
        let mut history_offsets = Vec::with_capacity(n_nodes);
        for _ in 0..n_nodes {
            let off = read_u64(&mut file)?;
            if off >= bytes {
                return Err(invalid("history offset beyond end of file"));
            }
            history_offsets.push(off);
        }
        Ok(Self {
            file,
            model,
            frame_offsets,
            history_offsets,
            flags,
            bytes,
        })
    }

    pub fn model(&self) -> &Mesh {
        &self.model
    }

    pub fn frame_count(&self) -> usize {
        self.frame_offsets.len()
    }

    pub fn file_len(&self) -> u64 {
        self.bytes
    }

    /// Field-presence union over all frames (file-level flags).
    pub fn has_temperatures(&self) -> bool {
        self.flags & FLAG_TEMPERATURES != 0
    }

    /// Field-presence union over all frames (file-level flags).
    pub fn has_velocities(&self) -> bool {
        self.flags & FLAG_VELOCITIES != 0
    }

    /// Read frame `i` (0-based) by seeking its recorded offset.
    pub fn frame(&self, i: usize) -> io::Result<OdbFrame> {
        let Some(&off) = self.frame_offsets.get(i) else {
            return Err(invalid(format!(
                "frame {i} out of range ({} frames)",
                self.frame_offsets.len()
            )));
        };
        let mut f = &self.file;
        f.seek(SeekFrom::Start(off))?;
        let flags = read_u32(&mut f)?;
        let time = read_f64(&mut f)?;
        let load_fraction = read_f64(&mut f)?;
        let n = read_u32(&mut f)? as usize;
        let displacements = read_triples(&mut f, n)?;
        let n = read_u32(&mut f)? as usize;
        let von_mises = read_f64s(&mut f, n)?;
        let n = read_u32(&mut f)? as usize;
        let nodal_von_mises = read_f64s(&mut f, n)?;
        let temperatures = if flags & FLAG_TEMPERATURES != 0 {
            let n = read_u32(&mut f)? as usize;
            Some(read_f64s(&mut f, n)?)
        } else {
            None
        };
        let velocities = if flags & FLAG_VELOCITIES != 0 {
            let n = read_u32(&mut f)? as usize;
            Some(read_triples(&mut f, n)?)
        } else {
            None
        };
        Ok(OdbFrame {
            time,
            load_fraction,
            displacements,
            von_mises,
            nodal_von_mises,
            temperatures,
            velocities,
        })
    }

    /// All frames in order (random access under the hood).
    pub fn frames(&self) -> impl Iterator<Item = io::Result<OdbFrame>> + '_ {
        let n = self.frame_offsets.len();
        (0..n).map(move |i| self.frame(i))
    }

    /// The |u| time series of one node (by model node index).
    pub fn history(&self, node_index: usize) -> io::Result<Vec<f64>> {
        let Some(&off) = self.history_offsets.get(node_index) else {
            return Err(invalid(format!(
                "history node {node_index} out of range ({} nodes)",
                self.history_offsets.len()
            )));
        };
        let mut f = &self.file;
        f.seek(SeekFrom::Start(off))?;
        read_f64s(&mut f, self.frame_offsets.len())
    }
}

fn read_model(f: &mut impl Read) -> io::Result<Mesh> {
    let n_nodes = read_u32(f)? as usize;
    let mut nodes = Vec::with_capacity(n_nodes);
    for _ in 0..n_nodes {
        let id = read_u32(f)?;
        nodes.push(solvix_core::Node {
            id,
            coords: [read_f64(f)?, read_f64(f)?, read_f64(f)?],
        });
    }
    let n_elements = read_u32(f)? as usize;
    let mut elements = Vec::with_capacity(n_elements);
    for _ in 0..n_elements {
        let id = read_u32(f)?;
        let kind = kind_from_code(read_u8(f)?)?;
        let n = read_u32(f)? as usize;
        let mut connectivity = Vec::with_capacity(n);
        for _ in 0..n {
            connectivity.push(read_u32(f)?);
        }
        elements.push(solvix_core::Element {
            id,
            kind,
            connectivity,
        });
    }
    Ok(Mesh { nodes, elements })
}

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::{Element, Node};

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "solvix_io_odb_{name}_{}.solvix",
            std::process::id()
        ))
    }

    fn demo_mesh() -> Mesh {
        Mesh {
            nodes: vec![
                Node {
                    id: 1,
                    coords: [0.0, 0.0, 0.0],
                },
                Node {
                    id: 2,
                    coords: [1.0, 0.5, -0.25],
                },
                Node {
                    id: 3,
                    coords: [1.0, 1.0, 0.0],
                },
                Node {
                    id: 4,
                    coords: [0.0, 1.0, 2.0],
                },
            ],
            elements: vec![
                Element {
                    id: 1,
                    kind: ElemKind::Cps4,
                    connectivity: vec![1, 2, 3, 4],
                },
                Element {
                    id: 2,
                    kind: ElemKind::Cps3,
                    connectivity: vec![2, 4, 3],
                },
            ],
        }
    }

    fn gnarly(i: usize) -> [f64; 3] {
        // Bit patterns round-trip tests care about: subnormals, -0.0,
        // extremes, PI.
        [
            f64::from_bits(0x0000_0000_0000_0001 + i as u64),
            -0.0,
            std::f64::consts::PI * i as f64 + f64::MIN,
        ]
    }

    fn demo_frame(i: usize, thermal: bool, vel: bool) -> OdbFrame {
        OdbFrame {
            time: i as f64 * 0.25,
            load_fraction: (i + 1) as f64 / 4.0,
            displacements: (0..4).map(|j| gnarly(i * 4 + j)).collect(),
            von_mises: vec![1e9 + i as f64, -2.5e-7],
            nodal_von_mises: vec![f64::MAX, 0.0, 1e-300, 42.0 + i as f64],
            temperatures: thermal.then(|| vec![100.0 - i as f64, 0.5, 0.0, -273.15]),
            velocities: vel.then(|| (0..4).map(|j| gnarly(i * 4 + j + 1)).collect()),
        }
    }

    fn write_frames(path: &Path, mesh: &Mesh, frames: &[OdbFrame]) -> io::Result<()> {
        let mut w = OdbWriter::create(path, mesh)?;
        for f in frames {
            w.append_frame(f)?;
        }
        w.finish()
    }

    #[test]
    fn round_trip_is_bit_exact() {
        let path = tmp("roundtrip");
        let mesh = demo_mesh();
        let frames = vec![
            demo_frame(0, true, false),
            demo_frame(1, true, false),
            demo_frame(2, true, true),
            demo_frame(3, true, true),
        ];
        write_frames(&path, &mesh, &frames).unwrap();
        let odb = Odb::open(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        assert_eq!(odb.model(), &mesh);
        assert_eq!(odb.frame_count(), 4);
        assert!(odb.has_temperatures());
        assert!(odb.has_velocities());
        for (i, got) in odb.frames().enumerate() {
            let got = got.unwrap();
            let want = &frames[i];
            assert_eq!(got.time.to_bits(), want.time.to_bits());
            assert_eq!(got.load_fraction.to_bits(), want.load_fraction.to_bits());
            assert_eq!(got, *want, "frame {i} must be bit-identical");
        }
    }

    #[test]
    fn appended_frames_read_back_in_order() {
        let path = tmp("append");
        let mesh = demo_mesh();
        let mut w = OdbWriter::create(&path, &mesh).unwrap();
        for i in 0..5 {
            w.append_frame(&demo_frame(i, false, false)).unwrap();
        }
        w.finish().unwrap();
        let odb = Odb::open(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(odb.frame_count(), 5);
        for (i, got) in odb.frames().enumerate() {
            assert_eq!(got.unwrap().time, i as f64 * 0.25);
        }
    }

    #[test]
    fn random_access_uses_correct_offsets() {
        let path = tmp("randaccess");
        let mesh = demo_mesh();
        let frames: Vec<OdbFrame> = (0..8).map(|i| demo_frame(i, false, false)).collect();
        write_frames(&path, &mesh, &frames).unwrap();
        let odb = Odb::open(&path).unwrap();

        // Offsets are strictly increasing and inside the file.
        let offs = &odb.frame_offsets;
        assert!(offs.windows(2).all(|w| w[0] < w[1]));
        assert!(offs.iter().all(|&o| o < odb.file_len()));

        // Read out of order: frame 5 first, then 0, then 7.
        assert_eq!(odb.frame(5).unwrap().time, 1.25);
        assert_eq!(odb.frame(0).unwrap().time, 0.0);
        assert_eq!(odb.frame(7).unwrap().load_fraction, 2.0);
        assert!(odb.frame(8).is_err());

        // History: node 3's |u| series matches manual computation.
        let want: Vec<f64> = frames
            .iter()
            .map(|f| {
                let d = f.displacements[3];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .collect();
        let got = odb.history(3).unwrap();
        for (g, w) in got.iter().zip(&want) {
            assert_eq!(g.to_bits(), w.to_bits());
        }
        assert!(odb.history(4).is_err());
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn bad_magic_and_version_rejected() {
        let path = tmp("badmagic");
        std::fs::write(&path, b"NOTSLVX010000000").unwrap();
        let err = Odb::open(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        let mut w = OdbWriter::create(&path, &demo_mesh()).unwrap();
        w.append_frame(&demo_frame(0, false, false)).unwrap();
        w.finish().unwrap();
        let mut raw = std::fs::read(&path).unwrap();
        // Version sits at bytes 8..12.
        raw[8..12].copy_from_slice(&99u32.to_le_bytes());
        std::fs::write(&path, &raw).unwrap();
        let err = Odb::open(&path).unwrap_err();
        assert!(err.to_string().contains("version 99"), "{err}");
        std::fs::remove_file(&path).unwrap();
    }

    /// 16x10 = 160 CPS4 quads, 187 nodes, 3 frames with every field.
    #[test]
    fn smoke_160_elements_quick() {
        let (nx, ny) = (16, 10);
        let nid = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;
        let mut mesh = Mesh::default();
        for j in 0..=ny {
            for i in 0..=nx {
                mesh.nodes.push(Node {
                    id: nid(i, j),
                    coords: [i as f64 / nx as f64, j as f64 / ny as f64, 0.0],
                });
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                mesh.elements.push(Element {
                    id: (j * nx + i + 1) as u32,
                    kind: ElemKind::Cps4,
                    connectivity: vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)],
                });
            }
        }
        assert_eq!(mesh.elements.len(), 160);
        let n = mesh.nodes.len();
        let frames: Vec<OdbFrame> = (0..3)
            .map(|i| OdbFrame {
                time: i as f64,
                load_fraction: (i + 1) as f64 / 3.0,
                displacements: vec![[0.001 * i as f64, -0.5, 0.0]; n],
                von_mises: vec![1e6 + i as f64; 160],
                nodal_von_mises: vec![2e6; n],
                temperatures: Some(vec![20.0 + i as f64; n]),
                velocities: Some(vec![[0.0, 1.5, 0.0]; n]),
            })
            .collect();
        let path = tmp("smoke");
        write_frames(&path, &mesh, &frames).unwrap();
        let odb = Odb::open(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(odb.frame_count(), 3);
        assert_eq!(odb.model().elements.len(), 160);
        let f2 = odb.frame(2).unwrap();
        assert_eq!(f2.von_mises.len(), 160);
        assert_eq!(f2.displacements.len(), n);
        assert_eq!(odb.history(0).unwrap().len(), 3);
    }
}
