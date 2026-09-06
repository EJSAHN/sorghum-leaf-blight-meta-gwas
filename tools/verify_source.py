"""Write or verify the source distribution's SHA-256 file inventory."""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path

EXCLUDED_DIRS = {'.git','__pycache__','.pytest_cache','.venv','venv','.mypy_cache','.ruff_cache'}
MANIFEST = 'SOURCE_MANIFEST_SHA256.txt'


def source_files(root: Path) -> dict[str, Path]:
    result = {}
    for path in root.rglob('*'):
        rel = path.relative_to(root)
        if any(x in EXCLUDED_DIRS for x in rel.parts) or not path.is_file():
            continue
        if path.name == MANIFEST or path.suffix in {'.pyc','.pyo'}:
            continue
        if path.is_symlink():
            raise ValueError('Source symlinks are not supported: '+rel.as_posix())
        result[rel.as_posix()] = path
    return result


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def write(root: Path) -> None:
    rows=['# SHA-256 source manifest; hashes cover the exact distributed file bytes.']
    rows += [digest(path)+'  '+name for name,path in sorted(source_files(root).items())]
    (root/MANIFEST).write_bytes(('\n'.join(rows)+'\n').encode('utf-8'))


def verify(root: Path) -> int:
    expected={}
    for line in (root/MANIFEST).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        sha, name=line.split(None,1)
        if name in expected or len(sha)!=64 or any(c not in '0123456789abcdefABCDEF' for c in sha):
            raise ValueError('Invalid source-manifest entry: '+name)
        expected[name]=sha.upper()
    files=source_files(root)
    if set(files)!=set(expected):
        raise ValueError('Source file-set mismatch: missing='+str(sorted(set(expected)-set(files)))+
                         '; extra='+str(sorted(set(files)-set(expected))))
    bad=[name for name,path in files.items() if digest(path)!=expected[name]]
    if bad:
        raise ValueError('Source checksum mismatch: '+', '.join(bad))
    return len(files)


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:
        write(args.root)
    print(f'Source manifest verified: {verify(args.root)} files.')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
