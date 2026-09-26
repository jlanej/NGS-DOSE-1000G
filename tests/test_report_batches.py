"""The sequencing batch of each generation of the complete trios (`report.trio_batches`, from NGS-PCA's release
batch): the 1000 Genomes 30x release sequenced its 698 related genomes, most of them trio children, after the original
2,504, which puts nearly all children in the later batch and most parents in the earlier one. Some parents were
sequenced with the children (F2 here); `shared` counts the families where the child's batch is also a parent's."""
from ngsdose.trios import Trio

from report.report import trio_batches


def test_the_batch_of_each_generation_is_counted():
    rows = [dict(sample=s, **{"ngspca.batch": b}) for s, b in (("C1", "698"), ("F1", "2504"), ("M1", "2504"), ("C2", "698"), ("F2", "698"), ("M2", "2504"))]
    b = trio_batches(rows, [Trio("C1", "F1", "M1"), Trio("C2", "F2", "M2"), Trio("C3", "F1", "M1")])
    assert b == dict(n=2, child={"698": 2}, parent={"2504": 3, "698": 1}, shared=1)
