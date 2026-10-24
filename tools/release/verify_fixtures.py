"""Verify the complete retained source-distribution TSV fixture set offline."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def verify():
    expected = json.loads((ROOT / 'release/fixture-sha256.json').read_text())
    folders = ['tests/fixtures', 'references/unrounded-probes',
               'references/legacy-oracles/paired', 'tie-order-requalification/fixtures']
    actual = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for folder in folders for p in (ROOT / folder).rglob('*.tsv')}
    # Research checkout contains additional paired oracle input data; the full
    # historical verifier is authoritative there. An archive must be exact.
    if (ROOT / '.git').exists():
        assert all(actual.get(name) == digest for name, digest in expected.items())
    else:
        assert actual == expected, 'Fixture membership or SHA-256 mismatch'
    print(f'PASS: {len(expected)} accepted TSV fixture SHA-256 identities')


if __name__ == '__main__':
    verify()
