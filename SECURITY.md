# Security Policy

## Supported versions

The maintainers investigate reports against the default branch and the most
recent tagged release. Older releases are reviewed as time and applicability
permit.

## Reporting a vulnerability

Please report possible security vulnerabilities through GitHub's private
security reporting:

<https://github.com/typezo/Solvix/security/advisories/new>

Do not open a public issue, pull request, or discussion for a suspected
vulnerability.

Include the information below where possible:

- description of the issue and its security impact;
- reproduction steps and a minimal input file;
- affected Solvix version or commit;
- environment details;
- any proposed mitigation.

The maintainers aim to acknowledge reports within three business days. Please
allow time for investigation and coordination before public disclosure. If you
do not receive a response through private reporting, follow up through the same
channel rather than opening a public issue.

## Scope and expectations

Solvix processes finite-element input files and executes solver workloads. The
maintainers therefore treat parser memory-safety violations, unauthorized file
access, excessive resource consumption from ordinary valid input, and unsafe
handling of untrusted data as security-relevant issues.

Ordinary crashes, numerical convergence problems, and unsupported solver
features should be reported as public bugs unless you can demonstrate a
security impact.
