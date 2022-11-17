import re
import json

def sequence_chunks(seq: str, peps, width: int):
    '''Show sequence with spans highlighting each mass spec peptide
    '''
    chunks = []
    pieces = []
    start = 0
    end = 0
    npeps = len(peps)
    if npeps == 0:
        pieces.append({'seq': seq[:width], 'is_pep': False})
        start = width
        end = width * 2
        chunks.append(pieces)
        while start < len(seq):
            pieces = [{'seq': seq[start:end], 'is_pep': False}]
            chunks.append(pieces)
            start = end
            end += width
        # print(json.dumps(chunks, indent = 4))
        return chunks
    loc = 0
    pep_idx = 0
    subseq = ''
    cur_pep = peps[0]
    end = cur_pep.location
    is_pep = False
    cutoff = False
    for ii, char in enumerate(seq):
        if ii % width == 0 and ii != 0:
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
                start = cur_pep.location
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
    # print(json.dumps(chunks, indent = 4))
    return chunks


def process_clustal_num(clustal: str, peptides, width: int):
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
    from .models import isoform_num
    sorted_acc_nums = sorted(seq_map.keys(), key = lambda x: 10000 if x == 'zzzz' else isoform_num(x))
    prots = []
    nchunks = 0
    for acc_num in sorted_acc_nums:
        seq = seq_map[acc_num]
        peps = (peptides
            .filter(prot = acc_num)
            .order_by('location')
        )
        chunks_this_seq = sequence_chunks(seq, peps, width)
        nchunks = max(len(chunks_this_seq), nchunks)
        prots.append({'acc_num': acc_num, 'chunks': chunks_this_seq})
    group_by_chunk = []
    old_chunk_ends = {}
    for ii in range(nchunks):
        curchunk = []
        for prot in prots:
            try:
                chunk = prot['chunks'][ii]
            except IndexError:
                chunk = [{'seq': ''}]
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
    # print(json.dumps(group_by_chunk, indent = 4))
    return {'header': header, 'chunks': group_by_chunk}