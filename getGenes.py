from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

args = argparse.ArgumentParser()
args.add_argument("--blastres", help="BLAST results file")
args.add_argument("--genome", help="Genome sequence file")
args.add_argument("--outfile", help="Output file")
args = args.parse_args()


atccfile = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
mappedSeqRec = []
print(atccfile)
blastRes = open(args.blastres, 'r', encoding="utf-8", newline="\n")
data = blastRes.readlines()
print(atccfile)
for identity in data:
    values = identity.split("\t")
    print(values)
    subSeq = values[1]
    seqStart = int(values[8])
    seqEnd = int(values[9])
    ident = float(values[2])
    print(ident)
    print(subSeq, seqStart, seqEnd)
    sequenceRec = atccfile[subSeq]
    print(sequenceRec)
    mappedSeq = ""
    sequence = str(sequenceRec.seq)
    if ident >= 90:
        if seqStart < seqEnd:
            mappedSeq = sequence[seqStart:seqEnd+1]
        else:
            mappedSeq = Seq(sequence[seqEnd:seqStart+1])
            mappedSeq = mappedSeq.reverse_complement()

    print(mappedSeq)
    if len(mappedSeq) != 0:
      mappedSeqRec.append(
        SeqRecord(id=values[0].split("_")[-1], seq=Seq(mappedSeq),
                  description=":".join([subSeq, str(seqStart), str(seqEnd)])))

SeqIO.write(mappedSeqRec, args.outfile, "fasta")
