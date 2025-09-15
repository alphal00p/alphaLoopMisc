#/usr/bin/env python3
import sys

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

start_reading = False
offset = 2
res = []
tot_time = 0.
tot_n_cuts = 0
n_graphs = 0
min_t = None
max_t = None
for line in open(sys.argv[1],'r').readlines():
    if not start_reading and not line.startswith('SG_ID SG_name n_cuts t_sample[ms]'):
        continue
    else:
        start_reading=True
    if offset > 0:
        offset -= 1
        continue
    if line.startswith('quit'):
        continue
    line = line.strip()
    #print([el for el in line.split(' ') if el != ''])
    sg_id, sg_name, n_cuts, t = [el for el in line.split(' ') if el != '']
    sg_id = int(sg_id)
    n_cuts = int(n_cuts)
    t = float(t)*1000
    if min_t is None or t < min_t:
        min_t = t
    if max_t is None or t > max_t:
        max_t = t
    tot_time += t
    tot_n_cuts += n_cuts
    n_graphs += 1
    res.append((sg_id, sg_name, n_cuts, t))


res = sorted(res, key= lambda el: el[-1])

template = r'\texttt{%-5s} & $%-2d$ & $%-5.1f$ '
for res_chunk in chunks(res, 6):
    for r in res_chunk[:-1]:
        print('{:<35s}'.format(template%(r[1],r[2],r[3]))+' | &')
    r = res_chunk[-1]
    print('{:<35s}'.format(template%(r[1],r[2],r[3])) + r'| \\')

print('')
print('')

avg_cuts = float(tot_n_cuts)/float(n_graphs)
avg_t = tot_time / float(n_graphs)
avg_t_per_cut = tot_time / float(tot_n_cuts)
print(r"""Avg. t: $\mathbf{%.1f}$\ \ \ 
Min. t: $\mathbf{%.1f}$\ \ \ 
Max. t: $\mathbf{%.1f}$\ \ \ 
Avg number of cuts: $\mathbf{%.1f}$\ \ \ 
Avg t per cut: $\mathbf{%.1f}$} \\"""%(
    avg_t,
    min_t,
    max_t,
    avg_cuts,
    avg_t_per_cut
))
print('')
print('')
print("Total time: %.1f"%tot_time)

