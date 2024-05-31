from libpymcr.utils import checkPath
import os
import subprocess

for v in ['R2021a', 'R2021b', 'R2022a', 'R2022b', 'R2023a', 'R2023b', 'R2024a']:
    print(f'Compiling for {v}')
    mlPath = checkPath(v)
    rv = subprocess.run([os.path.join(mlPath, 'bin', 'matlab'), '-batch', '"build_ctf; exit"'], capture_output=True)
    if rv.returncode != 0:
        print(rv.stdout.decode())
        print(rv.stderr.decode())
