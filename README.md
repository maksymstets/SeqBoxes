# SeqBoxes
Tool that highlights user-defined  subsequences (including protease/ restriction enzyme cleavage sites) of FASTA/FASTQ records in coloured boxes.

# seqboxes – User Guide


**Contents**

1. [What seqboxes does](#what-seqboxes-does)
2. [Installation and requirements](#installation-and-requirements)
3. [Quick start](#quick-start)
4. [Using the window](#using-the-window)
5. [Using the command line](#using-the-command-line)
6. [How sites are found, numbered and filtered](#how-sites-are-found-numbered-and-filtered)
7. [Results: the PNG files and the CSV file](#results-the-png-files-and-the-csv-file)
8. [Worked examples](#worked-examples)
9. [Limits and good to know](#limits-and-good-to-know)
10. [Troubleshooting](#troubleshooting)
11. [Running the tests and reading the code](#running-the-tests-and-reading-the-code)

## What seqboxes does

seqboxes turns FASTA or FASTQ sequences into pictures in which the sites chosen by user sit in coloured boxes, and it writes the box coordinates to a CSV file.

- **Protease cleavage sites.** Boxes the two residues on either side of every cleaved bond, using the 34 Expasy PeptideCutter rules that ship with Pyteomics.
- **Restriction-enzyme sites.** Boxes the recognition site of any Bio.Restriction enzyme in DNA, on either strand, in linear or circular sequences.
- **Your own subsequences.** Boxes every occurrence of any text you type, overlapping hits included.
- **Several things at once.** Each enzyme or subsequence has its own colour; boxes that overlap are nested so none is hidden.
- **Fewer boxes for broad enzymes.** Box only the first N sites, or only chosen site numbers. Numbering starts again in every sequence.
- **Multi-FASTA files.** The same settings apply to every sequence; each sequence is saved as its own PNG named after its SeqID.
- **A legend on every picture.** The enzyme names and their box colours are printed under the last row of the sequence.
- **Mixed protein and DNA files.** Proteases act on protein records and restriction enzymes on DNA records; the program tells them apart.

Start the program from a terminal and then work in a window. A command-line-only mode for batch runs also exists but is optional.

## Installation and requirements

Install three Python packages.

```bash
pip install biopython pyteomics matplotlib
```
The window needs tkinter, which Anaconda and the standard Windows and macOS Python installers already include. Without tkinter or a display only the command-line modes work.

Seqboxes was developed and tested with Python 3.12, Biopython 1.88 and Pyteomics 5.0.1. The number of restriction enzymes and the cleavage rules come from the installed packages, so other versions may list slightly different enzymes.

Two kinds of input file are read: FASTA (header line starts with `>`) and FASTQ (starts with `@`). Quality scores in FASTQ are ignored.

## Quick start

Run `python seqboxes.py` in a terminal, give it your FASTA file, and use the window that opens.

| What you type | What happens |
| --- | --- |
| `python seqboxes.py` | Asks for the file path in the terminal then opens the window. |
| `python seqboxes.py proteins.fasta` | Opens the window straight away with that file. |
| `python seqboxes.py proteins.fasta --enzyme trypsin --restriction EcoRI --circular` | Opens the window with these choices already filled in. |
| `python seqboxes.py proteins.fasta --enzyme trypsin --png-dir out --csv boxes.csv` | No window: writes the files and exits (batch mode). |

A first session in the window:

1. Start the program and give it the file.
2. Choose one sequence in the **Record(s)** list, or keep **(all records)**.
3. Click **+ Add protease** (or **+ Add restriction enzyme**), pick the enzyme and a colour.
4. Click **Render**. The picture appears on the right and a summary line appears under the controls.
5. Click **Save PNG(s)** and **Save CSV** and choose where the files should go.

Nothing is written to disk until saved by the user: closing the window leaves no files behind.

## Using the window

The controls are on the left and a scrollable preview is on the right; the picture changes only when clicked **Render**.

| Control | What it does |
| --- | --- |
| **Record(s)** | Chooses which sequence to show. With more than one sequence, **(all records)** comes first and draws every sequence one under another. |
| **Circular DNA** | Treats DNA records as circular: sites that cross the origin are found and drawn in two pieces. Protein records are not affected. A FASTA file does not say whether DNA is circular, so you decide here. |
| **Proteases** | **+ Add protease** adds a row. Choose the protease from the list of 34, click the colour swatch to pick a colour (or type a name such as `red` or a code such as `#1f77b4`), choose **box** (an outline) or **fill** (a highlight behind the letters), and use **x** to remove the row. |
| **Restriction enzymes** | **+ Add restriction enzyme** adds a row. The list shows commercially available enzymes, and you can type any enzyme name from the full database, for example `BsaI`; capital letters do not matter. An unknown name shows a warning. |
| **Which cleavage sites to box** | **All sites**, **First N sites** (the first N from the start of each sequence) or **Only site no.** (for example `3`, or `1,4,7-9`). One choice applies to every protease and restriction row. Invalid input shows a warning. |
| **Show site counts** | Opens a list with the number of sites for every sequence and enzyme. An entry such as `#3@45` means site number 3 starts at position 45. Look here before using **Only site no.** |
| **Subsequences to box** | Type any text, for example `LAG`. Every occurrence is boxed, overlapping ones too, whatever the letter case. **+ Add subsequence** adds another row. |
| **Skip sequences without matches…** | Appears when the file has several sequences. When ticked, sequences with no box are left out when PNG files are saved. |
| **Render** | Draws the picture and updates the summary line. |
| **Save PNG(s)** | With one sequence selected it asks for a file name. With several it asks for a folder and writes one PNG per sequence. |
| **Save CSV** | Asks for a file name and writes one CSV for all selected sequences. |

The summary line under the controls reports, for each row, how many sites or matches were found and in how many sequences, for example:

`3 sequence(s) shown. trypsin: 11 site(s) in 2/2 sequence(s); EcoRI: 1 site(s) in 1/1 sequence(s)`

If a row cannot act on the sequences shown, it says so (for example `EcoRI: skipped - no DNA sequences shown`). After saving, the same line shows where the file went.

The preview stacks all selected sequences in one picture. The saved PNG files are separate: each sequence gets its own picture and its own legend.

## Using the command line

The window opens unless you ask for a batch action; every option below either pre-fills the window or, for the batch actions, replaces it.

| Option | Meaning |
| --- | --- |
| `FILE` | FASTA or FASTQ file. If you leave it out, the program asks for it (not possible in batch mode). |
| `--record ID` | Use only this sequence, given as its SeqID or the start of it. Repeat the option for several. Default: all sequences. |
| `--enzyme NAME[=COLOUR[:box\|fill]]` | Box the cleavage sites of a protease in protein sequences. Repeat for several. |
| `--restriction NAME[=COLOUR[:box\|fill]]` | Box the recognition sites of a restriction enzyme in DNA sequences. Repeat for several. |
| `--rule SUBSEQUENCE=COLOUR[:box\|fill]` | Box every occurrence of a subsequence. Repeat for several. |
| `--max-sites N` | Box only the first N sites of each enzyme in each sequence. |
| `--site NUMBERS` | Box only these site numbers, for example `3` or `1,4,7-9`. Cannot be combined with `--max-sites`. |
| `--circular` | Treat DNA as circular. |
| `--count-sites` | Print the number of sites and their positions for each enzyme and sequence, then exit. |
| `--list-enzymes` | Print the 34 protease names and exit. |
| `--list-restriction` | Print all restriction enzyme names and exit. |
| `--png-dir DIR` | Batch mode: write one PNG per sequence into DIR (created if missing). |
| `--only-matches` | With `--png-dir`: skip sequences that have no box. |
| `--csv FILE` | Batch mode: write the CSV file. |
| `-h` | Show the built-in help. |

`--max-sites` and `--site` apply to every enzyme given and need at least one `--enzyme` or `--restriction`.

**Colours and styles.** A colour is any matplotlib colour name (`red`, `gold`) or a hex code (`#1f77b4`). The style is `box` (default) or `fill`. Enzymes given without a colour get these colours in order: red, blue, green, purple, orange, cyan, brown, pink.

**Enzyme names.** Capital letters do not matter and the start of a name is enough when it is unique, so `proteinase` means `proteinase k`. An ambiguous start such as `pepsin` is refused with the matching names listed. Put names that contain spaces in quotes.

**Notes.** If an enzyme cannot act on any selected sequence (a protease on DNA only, say), the program prints a `Note:` line and boxes nothing for it.

Examples:

```bash
# window, pre-filled
python seqboxes.py legumin.fasta --enzyme trypsin --enzyme "pepsin ph1.3=green"

# how many trypsin sites does each sequence have? (no window)
python seqboxes.py proteome.fasta --enzyme trypsin --count-sites

# only the second trypsin site of every sequence, pictures and CSV
python seqboxes.py proteome.fasta --enzyme trypsin --site 2 --png-dir out --csv boxes.csv

# restriction map of a plasmid, treated as circular
python seqboxes.py plasmid.fasta --restriction EcoRI --restriction BsaI=blue --circular --png-dir out
```

## How sites are found, numbered and filtered

Every row is a rule with its own colour, and each rule is applied to every sequence it can act on.

### Proteases

The 34 protease rules are the Expasy PeptideCutter specificities packaged in Pyteomics; trypsin also uses its published exceptions. Available names: arg-c, asp-n, bnps-skatole, caspase 1 to caspase 10, chymotrypsin high specificity, chymotrypsin low specificity, clostripain, cnbr, enterokinase, factor xa, formic acid, glutamyl endopeptidase, granzyme b, hydroxylamine, iodosobenzoic acid, lysc, ntcb, pepsin ph1.3, pepsin ph2.0, proline endopeptidase, proteinase k, staphylococcal peptidase i, thermolysin, thrombin, trypsin.

- **What is boxed.** The two residues either side of the cleaved bond, called P1 and P1'. The enzyme cuts between them.
- **Sequence ends.** A site needs both P1 and P1' to exist, so nothing is cut after the last residue. When a rule looks at a residue beyond a sequence end, that missing residue counts as allowed for "not X" conditions, so a site right at the start can still be found.
- **Thermolysin** follows the Expasy text, in which a proline at P2' blocks cleavage; the Expasy rule table leaves that out.
- **Not included.** LysN, TEV protease, neutrophil elastase and Asp-N with N-terminal Glu are not in Pyteomics' rule set.
- **Predictions only.** The rules summarise published specificities, and Expasy itself warns that experiments may differ.

### Restriction enzymes

Enzymes and recognition sequences come from Biopython's Bio.Restriction (REBASE data). Biopython 1.88 has 1,088 enzymes: 623 commercial and 465 not.

- **What is boxed.** The recognition site, found on either strand. Degenerate sites such as `GATNNNNATC` work.
- **Cut position.** The CSV column `cut_after` gives the number of the last base before the cut on the top strand.
- **Cut outside the sequence.** As in Bio.Restriction, a site in a linear sequence is not reported if its cut would fall outside the sequence.
- **No cut position known.** For 334 of the 1,088 enzymes REBASE lists only the recognition site. These are still boxed, and `cut_after` stays empty.
- **Two cuts per site.** 25 enzymes cut on both sides of their site; only the first cut is reported.
- **Circular DNA.** With **Circular DNA** ticked, a site that crosses the origin is found and drawn in two pieces, one at the end of the sequence and one at the start.

### Subsequences

A subsequence is matched literally, without wildcards, ignoring letter case. Overlapping hits are all boxed. Subsequences work on protein and DNA, and on circular DNA too.

### Numbering and limiting sites

The sites of one enzyme in one sequence are numbered 1, 2, 3 ... from the N-terminus (or the 5' end of DNA). Numbering starts again in every sequence and for every enzyme.

- **All sites** boxes everything.
- **First N sites** boxes sites 1 to N. A sequence with fewer sites gets all of them.
- **Only site no.** boxes the numbers you list. A sequence that has no such site gets no box.
- The number in the CSV column `site_number` is always the number among all sites, whatever the filter.

Use **Show site counts** in the window or `--count-sites` on the command line to see how many sites each sequence has before you choose.

### Protein or DNA

The program decides for every sequence on its own. A sequence is DNA if it contains only IUPAC DNA letters and at least 90 % of them are A, C, G, T or N; otherwise it is protein. Proteases act on protein, restriction enzymes on DNA, subsequences on both. A protein made only of such letters would be mistaken for DNA, but real proteins hardly ever are.

### Overlapping boxes

When several rules mark the same place, their boxes are drawn one inside the other with a small inset, so every colour stays visible.

## Results: the PNG files and the CSV file

The program produces one PNG per sequence and one CSV for the whole run, and it writes them only where you tell it to.

### Where the files go

| How you started | PNG files | CSV file |
| --- | --- | --- |
| Window, one sequence selected | **Save PNG(s)** asks for a file name; the file goes where you choose. | **Save CSV** asks for a file name. |
| Window, several sequences | **Save PNG(s)** asks for a folder and writes `<SeqID>.png` for each sequence in it. | **Save CSV** asks for a file name. |
| Command line with `--png-dir DIR` and `--csv FILE` | Written into DIR, which is created if it is missing. A relative path starts in the folder where you ran the command. | Written to FILE. |

After saving, the window's summary line shows the path, and batch mode prints a line such as `Wrote 3 PNG file(s) to out` and `Wrote boxes.csv (13 segments)`.

File names come from the SeqID with every character other than letters, digits, dot, underscore and hyphen replaced by an underscore: `gi|565165|gb|AAB31199.1|` is saved as `gi_565165_gb_AAB31199.1.png`. Two sequences with the same SeqID get `_2`, `_3` and so on, so nothing is overwritten by accident.

### The PNG

- **Title.** Only the SeqID, not the full FASTA header.
- **Rows.** 50 residues per row, in blocks of 10 separated by a space. The number above a block is the position of its last residue, so `10`, `20`, `30` and so on.
- **Boxes.** Outlines (or fills) in the colour of each rule; overlapping boxes are nested.
- **Legend.** Under the last row: a coloured box and the name of every protease, restriction enzyme and subsequence that acts on this sequence. A site filter is part of the name, for example `trypsin (first 3 sites)` or `EcoRI (site 2)`.
- **Grey notes** under the legend explain what a box means: P1 and P1' for proteases, recognition sites for restriction enzymes, and that sites may cross the origin in circular DNA.
- **Size.** At the default 300 dpi a picture is about 3,300 pixels wide and about 140 pixels tall per row of 50 residues. If a picture would become too large for a PNG file, the resolution is lowered automatically, to no less than 20 dpi.

### The CSV

The CSV has one row per drawn box piece. A site that wraps from one row of residues to the next, or across the origin of circular DNA, gives two rows. All positions start at 1 and include both ends.

| Column | Meaning |
| --- | --- |
| `record_id` | The SeqID of the sequence. |
| `enzyme` | Name of the protease or restriction enzyme; empty for a subsequence rule. |
| `enzyme_type` | `protease` or `restriction`; empty for a subsequence rule. |
| `site_number` | Number of the site among all sites of that enzyme in that sequence (numbering restarts in every sequence). Empty for subsequences. |
| `cut_after` | The cut falls after this position: the P1 residue of a protease, or the last base before the cut on the top strand of DNA. Empty for subsequences and for enzymes with no known cut position. |
| `subsequence` | The boxed letters, for example `KL` or `GAATTC`. |
| `colour`, `style` | Box colour as a hex code, and `box` or `fill`. |
| `match_start`, `match_end` | First and last position of the whole match. For a match across the origin `match_end` is smaller than `match_start`. |
| `segment_start`, `segment_end` | First and last position of the piece drawn in this row. |
| `row` | Row of the picture, counted from 1. |
| `x0_px`, `y0_px`, `x1_px`, `y1_px` | The box rectangle in pixels of that sequence's PNG: top-left corner (x0, y0) and bottom-right corner (x1, y1). The origin is the top-left corner of the image, x grows to the right and y downwards. |
| `topology` | `linear` or `circular` for DNA records; empty for proteins. |

The pixel columns match the PNG files written by the same program at the same resolution, so you can crop or annotate a box in other software. Two example rows (a trypsin site in a protein and an EcoRI site in DNA):

```csv
record_id,enzyme,enzyme_type,site_number,cut_after,subsequence,colour,style,match_start,match_end,segment_start,segment_end,row,x0_px,y0_px,x1_px,y1_px,topology
prot1,trypsin,protease,1,3,KL,#d62728,box,3,4,3,4,1,204.2,175.0,320.8,262.5,
mcs,EcoRI,restriction,1,1,GAATTC,#1f77b4,box,1,6,1,6,1,87.5,175.0,437.5,262.5,linear
```

## Worked examples

### Show only the first trypsin sites of many proteins

Broad enzymes put a box on a large part of a protein, so limit them. In the window: add **trypsin**, choose **First N sites**, enter `3`, click **Render**, then **Save PNG(s)** and pick a folder. Each protein gets a picture with at most three boxes, and the legend reads `trypsin (first 3 sites)`.

The same from the command line:

```bash
python seqboxes.py proteome.fasta --enzyme trypsin --max-sites 3 --png-dir out --csv boxes.csv
```

### Pick one particular site

First look up which numbers exist, then choose one.

```bash
python seqboxes.py proteome.fasta --enzyme trypsin --count-sites
python seqboxes.py proteome.fasta --enzyme trypsin --site 2 --png-dir out
```

The first command lists, for every sequence, entries such as `3, 24, 45` (the positions of sites 1, 2 and 3). The second boxes site 2 of every sequence; a sequence with fewer than two sites gets no box. In the CSV the `site_number` column then holds `2` on every row.

### A restriction map of a plasmid

Open a FASTA file with the plasmid, add **EcoRI** and **BamHI** with **+ Add restriction enzyme**, tick **Circular DNA** and click **Render**. If an EcoRI site spans the origin, its first letters at the end of the sequence and its last letters at the start are both boxed, and the grey note under the legend explains this. In the CSV that site has two rows with the same `match_start` (for example 119) and a smaller `match_end` (for example 3), and `topology` is `circular`.

### Protein and DNA in one file

A file may contain proteins and DNA together. Add `trypsin` and `EcoRI` and click **Render**: trypsin is applied to the protein records, EcoRI to the DNA records, and each picture's legend lists only the enzymes that act on it. The summary line says `skipped` for an enzyme that fits none of the sequences shown.

### An enzyme and your own motif

To see where trypsin cuts and where the motif `LAG` occurs, add `trypsin`, then type `LAG` in a subsequence row and set its style to **fill**. The legend lists both, and the motif appears as a highlight behind the letters while the cleavage boxes stay visible.

## Limits and good to know

- **Input formats.** Only FASTA and FASTQ are read. GenBank files and their circular/linear information are not.
- **Circular DNA is your choice.** A FASTA file does not say whether DNA is circular, so tick the box or use `--circular`. It applies to all DNA records in the run.
- **Protein or DNA is guessed.** See the section on how sites are found. A protein made only of A, C, G, T and N would count as DNA.
- **Sites, not peptides.** The picture shows where enzymes cut. It does not draw the peptides, missed cleavages or length limits.
- **Subsequences are literal.** No wildcards or patterns; use a restriction enzyme for degenerate DNA sites.
- **Restriction cut positions.** Enzymes with two cuts per site report only the first, and enzymes without a known cut have an empty `cut_after`.
- **Rules are predictions.** The protease rules summarise published specificities and may differ from an experiment.
- **Fixed layout.** Rows hold 50 residues in blocks of 10; this cannot be changed in the window.
- **Very long sequences.** Each residue is drawn individually, so sequences of hundreds of thousands of residues are slow and produce huge pictures. The resolution is lowered down to 20 dpi, and beyond about 300,000 residues in one sequence even that is not enough to keep the image within the PNG size limit.
- **Versions matter.** The list of restriction enzymes depends on your Biopython version and the protease rules on your Pyteomics version.
- **Overwriting.** Batch mode, and saving PNG files into a folder from the window, replace existing files of the same name without asking. Only the single-file save dialogs ask before replacing a file.

## Troubleshooting

Most problems are a missing package, a mistyped name or a rule that cannot act on the sequence type; the table lists the usual causes.

| Error | Cause and remedy |
| --- | --- |
| `ModuleNotFoundError: No module named 'pyteomics'` (or `Bio`, `matplotlib`) | A package is missing. Run `pip install biopython pyteomics matplotlib`. |
| `No module named 'tkinter'`, or no window opens | tkinter is missing or there is no screen (for example over SSH). Install `python3-tk` on Debian or Ubuntu, or use Anaconda. The batch options (`--png-dir`, `--csv`, `--count-sites`) work without a window. |
| `Error: Not a FASTA ('>') or FASTQ ('@') file.` | The first character of the file is neither `>` nor `@`. Check the file. |
| `Ambiguous protease 'pepsin'` | The start of the name fits several enzymes. Type more of it, for example `pepsin ph1.3`; the message lists the choices. |
| `Unknown restriction enzyme` or `Unknown protease` | Check the spelling. `--list-restriction` and `--list-enzymes` print the names your installation knows; a newer enzyme may be missing from an older Biopython. |
| A warning titled **Check input** | **First N** needs a whole number of at least 1, and **Only site no.** needs numbers such as `3` or `1,4,7-9`. The message says which. |
| `skipped - no DNA sequences shown` (or `protein`) | The enzyme does not fit the sequences shown: proteases need protein, restriction enzymes need DNA. |
| `0 site(s)` although you expected some | The site may truly be absent. Also check that **Only site no.** is not larger than the number of sites (**Show site counts**), that **Circular DNA** is ticked if the site crosses the origin, and remember that a restriction site is dropped in linear DNA when its cut would fall outside the sequence. |
| The picture did not change after you changed a setting | Click **Render**. |
| `Note: ... nothing will be boxed` on the command line | None of the selected sequences is of the type the enzyme needs. |
| `Record 'X' not found. Available: ...` | `--record` takes a SeqID or the start of one. The message lists the SeqIDs in the file. |


Every function and class in `seqboxes.py` has a docstring describing what it does, its arguments, what it returns and which errors it raises. Read them with `python -m pydoc seqboxes` or `help(seqboxes.render)` in a Python session.
