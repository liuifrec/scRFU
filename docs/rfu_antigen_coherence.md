# RFU Antigen-Evidence Coherence

`scrfu.tl.rfu_antigen_coherence()` asks whether exact, externally annotated
sequences within each RFU carry coherent antigen labels. It does not designate
an RFU as antigen-specific. The assignment policy is explicit: `nearest` uses
all assigned RFUs, while `threshold_pass` additionally requires `pass_thr`.

## Ambiguity and weighting

The default ambiguity policy is `fractional`. If one sequence has \(k\)
distinct antigen labels, it contributes \(1/k\) to each. Alternatives are
`exclude_ambiguous`, which omits sequences with more than one label, and
`multi_label`, which gives every observed label weight one. Evidence records
duplicating the same sequence-label pair do not increase antigen abundance.

`unique_sequence` weighting gives every distinct sequence one unit before
ambiguity allocation. `cell` weighting multiplies that contribution by the
number of selected receptor rows/cells carrying the sequence. Sequence and cell
match rates are both reported.

For RFU \(r\), let \(a_{rg}\) be the resulting abundance for antigen \(g\),
\(A_r=\sum_g a_{rg}\), and \(p_{rg}=a_{rg}/A_r\). Then:

- antigen richness is the number of labels with positive abundance;
- dominant fraction and antigen purity are
  \(\max_g(a_{rg})/A_r\);
- entropy is \(H_r=-\sum_g p_{rg}\log p_{rg}\);
- normalized entropy is \(H_r/\log K_r\) for richness \(K_r>1\), zero for one
  label, and undefined for no labels;
- ambiguous-sequence fraction is the fraction of labelled sequences with more
  than one distinct antigen label.

Coherence metrics are marked ineligible and returned as missing until
`min_matched_sequences` unique sequences contribute under the selected
ambiguity policy. VDJdb-matched, antigen-labelled, and coherence-contributing
sequence counts remain distinct. Counts, match rates, evidence records, scores,
and optional represented sample/donor counts remain available.

## Global descriptive measures

`global_antigen_coherence()` uses labelled unique sequences and reports:

- RFU purity and entropy averaged with matched-sequence weights;
- the fraction of labelled sequences in RFUs containing at least two distinct
  sequences;
- same-antigen pair fraction within RFUs, the mean label-overlap probability
  over all within-RFU labelled sequence pairs;
- between-RFU antigen sharing, the fraction of RFU pairs sharing at least one
  label;
- mutual information \(I(R;G)=\sum_{r,g}p(r,g)\log[p(r,g)/(p(r)p(g))]\);
- normalized mutual information \(I(R;G)/\sqrt{H(R)H(G)}\), when both marginal
  entropies are nonzero.

Fractional ambiguity weights define the joint distribution by default.
Undefined quantities remain missing rather than being silently replaced by
zero.
