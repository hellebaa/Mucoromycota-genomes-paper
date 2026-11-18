from pathlib import Path
import csv,json,hashlib
R=Path(__file__).resolve().parents[3]
V='hic_reference_20261118_v1';species=['gzMucCirc1','gzUmbVina1','ggMorZona1'];selected=set();sources=set()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 sources.add(p);return list(csv.DictReader(p.open(),delimiter='\t'))
def out(p,text):
 p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and p.read_text()!=text:raise RuntimeError('Refuse changing frozen release '+str(p))
 p.write_text(text);selected.add(p)
def table(p,rs):
 import io
 s=io.StringIO();w=csv.DictWriter(s,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs);out(p,s.getvalue())
assemblies=[];evidence=[];unresolved=[];tel=[]
for g in species:
 fasta=R/f'data/genomes/{g}/{g}.contigs.fasta';sizes={};c=None
 for l in fasta.open():
  if l.startswith('>'):c=l[1:].split()[0];assert c not in sizes;sizes[c]=0
  else:sizes[c]+=len(l.strip())
 assemblies.append(dict(genome_id=g,assembly_path=str(fasta.relative_to(R)),sha256=sha(fasta),contigs=len(sizes),bp=sum(sizes.values()),matrix_path=f'data/HiC/mcool/{g}.mcool',matrix_sha256=sha(R/f'data/HiC/mcool/{g}.mcool')))
 out(R/f'results/centromeres/tables/{V}__{g}__chrom.sizes',''.join(f'{c}\t{n}\n' for c,n in sizes.items()))
 rr=read(R/f'results/hic/tables/cross_species__harmonized_enrichment__{g}__refined_candidates.tsv');beds=[];targets=[];n=0
 for x in rr:
  c=x['chrom']
  if x['included'].lower()!='true':unresolved.append(dict(genome_id=g,contig=c,status='unresolved',reason='No retained candidate; potential assembly limitation remains a hypothesis'));continue
  n+=1;s=min(int(x['start']),int(x['broad_start']));e=max(int(x['end']),int(x['broad_end']));assert 0<=s<e<=sizes[c]
  name=f'{g}__{c}__candidate';beds.append(f'{c}\t{s}\t{e}\t{name}\t0\t.\n');targets.append(f"{c}\t{x['start']}\t{x['end']}\t{name}__enrichment_target\t0\t.\n")
  evidence.append(dict(genome_id=g,contig=c,candidate=name,broad_start0=s,broad_end=e,enrichment_start0=x['start'],enrichment_end=x['end'],anchor0=x['center'],status='curated_putative',selection=x.get('selection',''),contact_evidence='Hi-C inspected; same-data enrichment/refinement, not independent validation',prior_evidence='Published Baalsrud et al. G3 BLAST-based annotation using MucCirc3 CEN sequence' if g=='gzMucCirc1' else 'No published sequence anchor used',sequence_evidence='See shortlisted motif BED and source audit' if g=='gzMucCirc1' else 'LINE-associated motif evidence is heterogeneous; not required for inclusion' if g=='gzUmbVina1' else 'Repeat-associated motifs heterogeneous; no universal sequence signature established',notes=x.get('curation_note',x.get('note','')),source=str((R/f'results/hic/tables/cross_species__harmonized_enrichment__{g}__refined_candidates.tsv').relative_to(R))))
 assert n=={'gzMucCirc1':10,'gzUmbVina1':12,'ggMorZona1':20}[g]
 if g=='ggMorZona1':unresolved.append(dict(genome_id=g,contig='contig_18',status='excluded',reason='User excluded uncertain assembly/contact pattern'))
 if g=='gzUmbVina1':
  for c in ['contig_13','contig_14','contig_15','contig_16']:unresolved.append(dict(genome_id=g,contig=c,status='no_retained_candidate',reason='No clear candidate retained in curation'))
 out(R/f'results/centromeres/tables/{V}__{g}__candidates.bed',''.join(beds));out(R/f'results/centromeres/tables/{V}__{g}__enrichment_targets.bed',''.join(targets))
 mort=g=='ggMorZona1';stem='mortierellales_T6AGGG_20261117_retry2' if mort else 'telomere'
 bp=R/('results/assembly_structure/tables/'+(stem+'__terminal_hits.bed' if mort else 'telomere_hits_multi_motif.bed'));sources.add(bp)
 calls=[];strict=[]
 for l in bp.open():
  v=l.strip().split('\t')
  if len(v)<6 or v[3].split('|')[0]!=g:continue
  c=v[0];s=int(v[1]);e=int(v[2]);assert 0<=s<e<=sizes[c];side=v[3].split('|')[1];dist=s if side.startswith('left') else sizes[c]-e;copies=int(v[4]);orientation=v[5]==('-' if side.startswith('left') else '+');keep=orientation and dist<=100 and copies>=5
  calls.append('\t'.join(v)+'\n')
  if keep:strict.append('\t'.join(v)+'\n')
  tel.append(dict(genome_id=g,contig=c,start0=s,end=e,name=v[3],strand=v[5],copies=copies,distance_to_end_bp=dist,expected_orientation=orientation,stringent_terminal_call=keep,source=str(bp.relative_to(R))))
 out(R/f'results/assembly_structure/tables/{V}__{g}__telomeres_all.bed',''.join(calls));out(R/f'results/assembly_structure/tables/{V}__{g}__telomeres_stringent.bed',''.join(strict))
# Molecular motif positions remain separate from broad interval bounds.
mp=R/'results/centromeres/tables/gzMucCirc1__motif_position_audit_20261030.tsv';mr=read(mp)
chosen=[x for x in mr if x['pooled_AT_shortlist']=='True' and x['retained_reference_on_contig']=='True'];assert len(chosen)==10
out(R/f'results/centromeres/tables/{V}__gzMucCirc1__shortlisted_motifs.bed',''.join(f"{x['chrom']}\t{x['start']}\t{x['end']}\t{x['chrom']}__shortlisted_CEN_motif\t0\t{x['strand']}\n" for x in chosen))
table(R/f'results/centromeres/tables/{V}__evidence.tsv',evidence);table(R/f'results/centromeres/tables/{V}__unresolved.tsv',unresolved);table(R/f'results/assembly_structure/tables/{V}__telomere_evidence.tsv',tel);table(R/f'results/hic/tables/{V}__assemblies.tsv',assemblies)
# Only compact exports are released; raw analytical inputs stay local for now.
readme=R/'analysis/hic/README.md'
out(R/f'results/hic/logs/{V}__README.md',readme.read_text())
provenance={'release':V,'candidate_counts':{'gzMucCirc1':10,'gzUmbVina1':12,'ggMorZona1':20},'builder_sha256':sha(Path(__file__)),'source_sha256':{str(p.relative_to(R)):sha(p) for p in sorted(sources)},'coordinate_convention':'0-based half-open','matrix_status':'original; filtered/rebalanced analysis pending'}
out(R/f'results/hic/logs/{V}__provenance.json',json.dumps(provenance,indent=2,sort_keys=True)+'\n')
checksum=R/f'results/hic/logs/{V}__checksums.tsv'
# A checksum manifest cannot include its own hash; the Git manifest covers it.
checksum_rows=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(selected)]
table(checksum,checksum_rows)
assert len(selected)==23, len(selected)
manifest=R/'analysis/hic/config/reference_v1_files.tsv';manifest.parent.mkdir(parents=True,exist_ok=True)
with manifest.open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','bytes','sha256'])
 for p in sorted(selected):w.writerow([str(p.relative_to(R)),p.stat().st_size,sha(p)])
with (manifest.parent/'reference_v1_sources.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','bytes','sha256','availability'])
 for p in sorted(sources):w.writerow([str(p.relative_to(R)),p.stat().st_size,sha(p),'local; not part of compact NeLS release'])
print('Compact release:',len(selected),'files;',sum(p.stat().st_size for p in selected),'bytes')
