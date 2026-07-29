# ChatGPT prompt for reproducing the Hammock workflow diagram

Use the following prompt to recreate the current version of the paper workflow figure. The target is a polished, publication-quality scientific infographic rather than a decorative illustration.

## Image-generation prompt

Create a high-resolution scientific workflow figure for a genomics methods paper. Use a clean, flat, vector-infographic style on a white background, suitable for a journal manuscript. The final image should be landscape orientation with an approximately 16:9 aspect ratio. Use crisp sans-serif typography, restrained color, rounded rectangular panel borders, consistent arrows and icons, and enough whitespace that the figure remains legible when reduced for publication.

The figure is a single multipanel figure with three labeled panels: **A**, **B**, and **C**. Panel A is a narrow vertical panel on the left. Panels B and C are two wide horizontal workflow lanes stacked on the right, with B above C. The complete composition should be centered and scaled to fill the canvas evenly; do not leave a large empty region on the right. Beneath the panels, include one centered takeaway banner spanning most of the figure width.

At the very top, centered, place the bold title:

**Figure 2. Hammock workflow**

Directly underneath, in smaller italic text, place the subtitle:

**Two complementary sketch-based representations of genomic interval sets**

Do not include a separate comparison table on the right side.

---

## Panel A: Inputs and analysis barriers

Create a tall rounded rectangular panel with a thin dark gray or black border and a very light gray background. Label the panel at the top in bold:

**A. Inputs and analysis barriers**

Below it, add the subheading:

**Example BED files across references**

Show three example BED-file rows. Each row should contain:

- a white document icon with a folded corner and a blue label reading `BED`;
- the filename in black;
- the reference assembly in blue;
- the label `chr1`;
- a miniature chromosome or genome-track schematic with alternating gray, black, and white segments and a few blue interval blocks;
- an ellipsis at the right edge of the track.

Use these exact examples:

1. `Sample_A.bed` — `hg19`
2. `Sample_B.bed` — `hg38`
3. `Sample_C.bed` — `CHM13`

Separate the next section with a light dashed horizontal divider.

### A1. All-pairs scaling

Add a dark circular number marker containing `1`, followed by the bold label:

**All-pairs scaling**

Under it, show the formula:

`N(N−1)/2`

Place a small dense network diagram beside the formula. The network should show several small BED-file icons connected to one another by many lines, visually representing all pairwise comparisons. Add an ellipsis to imply that the number of files can be much larger.

Separate the next section with another light dashed divider.

### A2. Reference-specific coordinates

Add a dark circular number marker containing `2`, followed by the bold label:

**Reference-specific coordinates**

Show two chromosome-track schematics aligned vertically:

- the upper track labeled `hg19 chr1`;
- the lower track labeled `hg38 chr1`.

Place a vertical broken or dashed connector between the tracks with a prominent red `X` to show that direct coordinate overlap is not valid across references. Keep this schematic simple and immediately understandable.

At the bottom of Panel A, add a small rounded gray callout box. Inside it, place a simple dark gray database-stack, archive, or repository icon—not a people icon. Next to the icon, use the exact two-line text:

**Public interval collections are large**  
**and distributed across references.**

The icon should be compact and secondary to the text.

From the right edge of Panel A, draw two thick arrows that split into the two workflows:

- a blue arrow pointing into Panel B;
- a green arrow pointing into Panel C.

These arrows should clearly show that interval mode and sequence mode are parallel representations of the same BED-file inputs, not sequential steps.

---

## Panel B: Interval mode

Create a wide rounded rectangular panel with a thin blue border and an extremely light blue background. Label it in bold blue:

**B. Interval mode (within-reference coordinate similarity)**

Arrange four numbered steps from left to right. Use blue circular step markers and black arrows between steps.

### B1. BED intervals

Heading:

**BED intervals**

Under it, show:

- the label `chr1 (hg19)`;
- a small chromosome schematic;
- several horizontal genomic tracks with blue interval blocks at different positions;
- an ellipsis indicating additional intervals or files.

This should look like ordinary BED interval annotations on a chromosome.

### B2. Covered genomic positions

Heading:

**Covered genomic positions**

Under it, again show `chr1 (hg19)`, but now visualize the intervals as the union of covered base-pair positions. Use a dense row of narrow blue marks or tiles spanning the covered regions. Beneath the track, show a compact set notation example:

`{101,102,103,...}`

Below that, place the blue explanatory text:

**Elements sketched:**  
**covered positions**

The visual should make clear that the sketch represents covered genomic positions rather than treating interval records as indivisible tokens.

### B3. HLL sketch

Heading:

**HLL sketch**

Under it, place the label:

**Reusable**  
**fixed-size sketch**

Show a compact HyperLogLog register-array icon: a long horizontal series of small outlined cells containing varying blue fill heights or intensities. Label a few positions above the array as:

`1   2   3   ...   m−1   m`

Below the array, add the blue bullet points:

- **built once**
- **one sketch per BED file**

The emphasis should be that each input BED file is read and summarized once, after which the fixed-size sketch can be reused for many comparisons.

### B4. Similarity matrix / all-pairs comparison

Heading:

**Similarity matrix /**  
**all-pairs comparison**

Show several small blue sketch-array icons stacked vertically on the left, with an ellipsis indicating many sketches. Draw arrows from those sketches into a compact square similarity heatmap. The matrix should use pale-to-dark blue cells with a dark diagonal and lighter off-diagonal values.

Near the matrix, add the blue italic text:

**Fast all-pairs comparison**  
**within a shared reference**

At the bottom of Panel B, include a narrow light-blue footer strip spanning the panel. Center the following text within it, separated by a vertical divider:

**Question answered: Do these files cover similar genomic locations?**  |  **Requires a common coordinate system.**

---

## Panel C: Sequence mode

Create a wide rounded rectangular panel directly below Panel B. Use a thin green border and an extremely light green background. Label it in bold green:

**C. Sequence mode (cross-reference sequence similarity)**

Arrange five numbered steps from left to right. Use green circular step markers and black arrows between steps.

### C1. BED + native reference FASTA

Heading:

**BED + native**  
**reference FASTA**

Show a blue BED-file icon with a chromosome track containing green highlighted intervals. Below or beside it, place a plus sign and a white FASTA document icon with a green label reading `FASTA`.

Next to the FASTA icon, list:

- `hg19 FASTA`
- `hg38 FASTA`
- `CHM13 FASTA`

The visual must clearly communicate that each BED file is paired with the reference FASTA corresponding to the coordinate system in which that BED file was defined.

### C2. Extract interval sequences

Heading:

**Extract interval**  
**sequences**

Show a chromosome track with one interval highlighted in green. Use dashed guide lines from that interval to a short nucleotide string such as:

`ACCTGAGT...TCCGA`

Below it, show a small FASTA-like file labeled:

`intervals.fa`

Inside the file, include a few example entries:

`>interval_1`  
`ACCTGAGT...TCCGA`

`>interval_2`  
`TGGACCTA...AGGCT`

`>interval_3`  
`GCTTACGA...TCGAA`

`...`

Make it visually clear that one BED file becomes one collection of interval-derived sequences, not one separate sketch for every interval.

### C3. Minimizers

Heading:

**Minimizers**

Show a nucleotide sequence such as:

`ACCTGAGTCGATCCGA`

Above or below the sequence, draw several overlapping dashed window boxes to represent sliding windows. Within those windows, highlight one representative k-mer per window. Draw arrows from selected windows down to compact green-outlined k-mer boxes labeled, for example:

- `ACCTG`
- `GTCGA`
- `TCCGA`
- `...`

Below these, add the green text:

**Representative k-mers**  
**selected by minimizers**

Do not include parameter values, flanking-k-mer options, command-line syntax, or implementation details beyond the conceptual minimizer selection.

### C4. HLL sketch

Heading:

**HLL sketch**

Under it, use the same sketch icon design as in Panel B, but with green fills. Add the label:

**Reusable**  
**fixed-size sketch**

Label register positions above the array as:

`1   2   3   ...   m−1   m`

Below the sketch, add the green bullet points:

- **built once**
- **one sketch per BED file**

This repeated visual style should emphasize that both modes share the same sketch-based compression framework even though they sketch different elements.

### C5. Cross-reference similarity / clustering

Heading:

**Cross-reference similarity /**  
**clustering**

Show nine small green sketch-array icons arranged in three groups, with the following labels:

- `Heart (hg19)`
- `Heart (hg38)`
- `Heart (CHM13)`
- `Liver (hg19)`
- `Liver (hg38)`
- `Liver (CHM13)`
- `Lung (hg19)`
- `Lung (hg38)`
- `Lung (CHM13)`

To the right, connect them with a compact dendrogram or clustering tree. The tree should group samples by tissue rather than reference:

- all Heart samples in a red-labeled `Heart cluster`;
- all Liver samples in a blue-labeled `Liver cluster`;
- all Lung samples in a purple-labeled `Lung cluster`.

Keep the dendrogram illustrative rather than data-heavy. Its role is to show the intended cross-reference biological comparison, not to present measured values.

At the bottom of Panel C, include a narrow light-green footer strip spanning the panel. Center the following text within it, separated by a vertical divider:

**Question answered: Do these files contain similar underlying sequence?**  |  **No direct coordinate overlap required.**

---

## Bottom takeaway banner

Below Panels A, B, and C, add one wide rounded rectangular callout banner centered horizontally. It should span most of the figure width and have a thin muted orange or gold border with a very pale cream background.

Do **not** place a star or any other icon in this banner.

Center the takeaway text both horizontally and vertically. Use one readable line of bold text, sized to fill the banner without crowding:

**One tool, two complementary similarity spaces: coordinates within references and sequence across references.**

Color only the key phrases:

- `coordinates within references` in blue;
- `sequence across references` in green;
- all other text in black.

The text should be visually centered and evenly balanced within the callout box.

---

## Composition and styling constraints

- Treat the graphic as one figure with three labeled panels: A, B, and C.
- Use a clean journal-figure aesthetic rather than a marketing infographic.
- Keep the entire composition centered and scaled to fill the 16:9 canvas evenly.
- Minimize unused white space, especially on the right side.
- Preserve generous internal whitespace between workflow steps.
- Use blue consistently for interval mode and green consistently for sequence mode.
- Use neutral black, gray, and white for shared inputs and arrows.
- Use red only for the coordinate mismatch `X` and the Heart cluster label.
- Use purple only for the Lung cluster label.
- Maintain consistent icon sizes and line weights across panels.
- Keep text large enough to remain legible after reduction in a manuscript.
- Make all arrows directional and unambiguous.
- Do not include the previously proposed right-side mode-comparison table.
- Do not show sequence mode as downstream of interval mode; both must branch independently from Panel A.
- Do not describe sequence mode as reference-free. It uses each BED file's native reference FASTA but enables comparison without shared coordinates.
- Do not label sequence-mode output as cross-reference overlap. Use `cross-reference similarity` or `sequence-content similarity`.
- Do not imply that sequence similarity establishes coordinate homology.
- Do not add numerical results, benchmark values, software commands, or parameter defaults.
- Avoid decorative illustrations, gradients beyond subtle heatmap shading, 3D effects, drop shadows, or photographic elements.

The final image should read as a conceptual map of the paper: Panel A introduces the two barriers, Panel B shows scalable coordinate-based comparison within one reference, and Panel C shows sequence-content comparison across references.