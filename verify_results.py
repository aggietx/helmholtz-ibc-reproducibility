"""Check the saved manuscript rows without running a PDE solve."""
from pathlib import Path
import re

root = Path(__file__).resolve().parent / 'results'
pattern = re.compile(r'TRACE_RESULT N=(\d+) ell=(\d+) method=(\w+) PC=(\d+) flag=(\d+) residual=([\deE.+-]+)')
grids = (64,128,256,512)
for filename, levels, expected in (
    ('table1_asymptotic.log', (2,2,3,3), [(3,3)]*4),
    ('table2_hthreehalves.log', (3,3,3,4), [(2,2)]*4),
    ('table3_ppw20_2d.log', (1,1,1,1), [(4,5),(5,5),(5,5),(5,6)])):
    records = {(int(n),int(ell),m):(int(pc),int(flag),float(r))
               for n,ell,m,pc,flag,r in pattern.findall((root/filename).read_text())}
    for n,ell,counts in zip(grids,levels,expected):
        for method,count in zip(('fgmres20','stationary'),counts):
            pc,flag,res = records[n,ell,method]
            assert pc == count and flag == 0 and res <= 1e-6
    print(filename, 'PASS')
for n,count in zip(grids,(4,4,5,5)):
    rows = re.findall(r'OUTER_RESULT (\w+) N=(\d+) PC=(\d+) reason=(-?\d+) residual=([\deE.+-]+)',
                      (root/f'table4_ppw20_3d_{n}.log').read_text())
    assert len(rows) == 2
    for method,mesh,pc,reason,res in rows:
        assert int(mesh)==n and int(pc)==count and int(reason)>0 and float(res)<=1e-6
    print(f'table4 N={n} PASS')
