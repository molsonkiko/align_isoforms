import re
# from gorp import gprint

def sequence_chunks(seq: str, peptides):
    '''Show sequence with spans highlighting each mass spec peptide
    '''
    chunks = []
    start = 0
    end = 0
    for pep in peptides:
        end = pep.location
        chunks.append({'seq': seq[start:end], 'is_pep': False})
        chunks.append({'seq': pep.peptide, 'is_pep': True, 'loc': end})
        start = end
    chunks.append({'seq': seq[end:], 'is_pep': False})
    return chunks


def process_clustal_num(clustal: str, peptides):
    chunks = re.split('\n{2,3}', clustal)
    header = chunks[0]
    seq_map = {}
    stars = ''
    for chunk in chunks[1:]:
        lines = chunk.split('\n')
        for line in lines:
            if not line:
                continue
            if line[0] == ' ':
                stars += line[14:]
                continue
            acc_num, seq = line.split()[:2]
            seq_map.setdefault(acc_num, '')
            seq_map[acc_num] += seq
    seq_map['zzzz'] = stars
    sorted_acc_nums = sorted(seq_map.keys())
    prots = []
    nchunks = 0
    for acc_num in sorted_acc_nums:
        seq = seq_map[acc_num]
        peps = (peptides
            .filter(prot = acc_num)
            .order_by('location')
        )
        npeps = len(peps)
        chunks = []
        pieces = []
        if npeps == 0:
            pieces.append({'seq': seq[:60], 'is_pep': False})
            start = 60
            end = 120
            while end < len(seq):
                chunks.append(pieces)
                pieces = [{'seq': seq[start:end], 'is_pep': False}]
                start = end
                end += 60
            chunks.append([{'seq': seq[end:], 'is_pep': False}])
            prots.append({
                'acc_num': acc_num,
                'chunks': chunks
            })
            nchunks = max(nchunks, len(chunks))
            continue
        loc = 0
        pep_idx = 0
        subseq = ''
        cur_pep = peps[0]
        start = 0
        end = cur_pep.location
        is_pep = False
        cutoff = False
        for ii, char in enumerate(seq):
            if ii % 60 == 0 and ii != 0:
                pep = {
                    'seq': subseq, 
                    'is_pep': is_pep,
                }
                if is_pep:
                    pep['loc'] = start
                    pep['pep_num'] = 0
                cutoff = True
                pieces.append(pep)
                subseq = ''
                chunks.append(pieces)
                pieces = []
            if char == '-':
                subseq += char
                continue
            if loc >= end:
                if loc == end:
                    pep = {
                        'seq': subseq, 
                        'is_pep': is_pep,
                    }
                    if is_pep:
                        pep['loc'] = start
                        pep['pep_num'] = int(cutoff)
                    cutoff = False
                    pieces.append(pep)
                    subseq = ''
                    start = end
                else:
                    start = loc - 1
                if is_pep:
                    pep_idx += 1
                    if pep_idx < npeps:
                        cur_pep = peps[pep_idx]
                        end = cur_pep.location
                    else:
                        end = len(seq)
                else:
                    end = cur_pep.location + len(cur_pep.peptide)
                is_pep = not is_pep
            loc += 1
            subseq += char
        pep = {
            'seq': subseq, 
            'is_pep': is_pep,
        }
        if is_pep:
            pep['loc'] = start
            pep['pep_num'] = int(cutoff)
        pieces.append(pep)
        chunks.append(pieces)
        nchunks = max(nchunks, len(chunks))
        prots.append({
            'acc_num': acc_num,
            'chunks': chunks,
        })
    group_by_chunk = []
    old_chunk_ends = {}
    for ii in range(nchunks):
        curchunk = []
        for prot in prots:
            try:
                chunk = prot['chunks'][ii]
            except IndexError:
                continue
            acc_num = '' if prot['acc_num'] == 'zzzz' else prot['acc_num']
            old_chunk_ends.setdefault(acc_num, 0)
            old_chunk_end = old_chunk_ends[acc_num]
            chunk_end = old_chunk_end + sum(len(p['seq']) - p['seq'].count('-') for p in chunk)
            is_prot_chunk = len(acc_num) > 0
            curchunk.append({
                'acc_num': acc_num,
                'chunk': chunk,
                'chunk_end': chunk_end,
                'is_prot_chunk': is_prot_chunk,
            })
            old_chunk_ends[acc_num] = chunk_end
        group_by_chunk.append(curchunk)
    # gprint.gprint(group_by_chunk)
    return {'header': header, 'chunks': group_by_chunk}