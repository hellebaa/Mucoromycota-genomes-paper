#!/usr/bin/env python3
"""Fetch a checked repository-relative NeLS manifest in one rsync session."""
import csv,hashlib,subprocess,sys,tempfile
from pathlib import Path,PurePosixPath
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader(Path(sys.argv[1]).open(),delimiter='\t'))
for x in rows:
 p=PurePosixPath(x['path'])
 if p.is_absolute() or '..' in p.parts or p.parts[0] not in ('data','results') or '\n' in x['path']:raise ValueError('Unsafe manifest path')
 if len(x['sha256'])!=64:raise ValueError('Invalid checksum')
with tempfile.NamedTemporaryFile(mode='w') as f:
 f.write(''.join(x['path']+'\n' for x in rows));f.flush()
 subprocess.run(['rsync','-av','--files-from='+f.name,'nels.storage:Projects/NMBU_threeD_Sandve_2024/Mucor/',str(root)+'/'],check=True)
for x in rows:
 p=root/x['path'];h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 if h.hexdigest()!=x['sha256'] or p.stat().st_size!=int(x['bytes']):raise RuntimeError('Checksum/size mismatch: '+x['path'])
print('Verified',len(rows),'downloaded files')
