# DNA Compression Research

Analyzing different compression algorithms for genomic sequencing data

The  Jupyter Notebooks in the repository contain further information

- **`dna-notes.ipynb`** - Notes made from different research papers and websites
- **`dna-sample-analysis.ipynb`** - Analyzes the protein sequence in `data/sample.fa` and gives information about the frequency of different nucleobases in the sample FASTA sequence
- **`dynamic-markov-compression.ipynb`** - Implements ***first-order dynamic Markov compression*** algorithm proposed in the paper by G. V. Cormack and R. N. Horspool [<sup>[1]</sup>](#1) to find the compression ratio on the genomic data of a dengue virus `data/dengue.fa`

## Algorithms

### Dynamic Markov Compression

Dynamic Markov Compression (DMC) uses predictive arithmetic coding similar to the prediction by partial matching (PPM), except that the input is predicted one bit at a time (rather than one byte at a time). DMC has a good compression ratio and moderate speed, similar to PPM, but requires somewhat more memory and is not widely implemented.

Implementing ***Markov modeling*** and ***Guazzo's arithmetic coding*** together provides a powerful data compression method. Its advantage is that it is adaptive, i.e. messages can be encoded and decoded with a single pass through the data.

All compression algorithms rely on a priori assumption of the data. The new direction taken here is that an algorithmic attempt is used to discover a Markov chain model that describes the data. If it can be created using the first part of the data then it can be used to predict forthcoming characters.

The decoder will generate the same statistical properties that the encoder is using for the data as the decoder is decoding it.

## Getting started

- *(Optional)* Install [Anaconda](https://www.anaconda.com/distribution/#download-section)
- Install [JupyterLab](https://jupyter.org/install)
- Clone the repository
- Create a virtual environment using `requirements/requirements.txt` or `requirements/conda-requirements.txt`
- Open the notebooks using JupyterLab

## Results 

#### Contents of the sample file 

`head -n 10 data/dengue.fasta`
```
>NC_001477.1 Dengue virus 1, complete genome
AGTTGTTAGTCTACGTGGACCGACAAGAACAGTTTCGAATCGGAAGCTTGCTTAACGTAGTTCTAACAGT
TTTTTATTAGAGAGCAGATCTCTGATGAACAACCAACGGAAAAAGACGGGTCGACCGTCTTTCAATATGC
TGAAACGCGCGAGAAACCGCGTGTCAACTGTTTCACAGTTGGCGAAGAGATTCTCAAAAGGATTGCTTTC
AGGCCAAGGACCCATGAAATTGGTGATGGCTTTTATAGCATTCCTAAGATTTCTAGCCATACCTCCAACA
GCAGGAATTTTGGCTAGATGGGGCTCATTCAAGAAGAATGGAGCGATCAAAGTGTTACGGGGTTTCAAGA
AAGAAATCTCAAACATGTTGAACATAATGAACAGGAGGAAAAGATCTGTGACCATGCTCCTCATGCTGCT
GCCCACAGCCCTGGCGTTCCATCTGACCACCCGAGGGGGAGAGCCGCACATGATAGTTAGCAAGCAGGAA
AGAGGAAAATCACTTTTGTTTAAGACCTCTGCAGGTGTCAACATGTGCACCCTTATTGCAATGGATTTGG
GAGAGTTATGTGAGGACACAATGACCTACAAATGCCCCCGGATCACTGAGACGGAACCAGATGACGTTGA
```

#### Markov chain model created using the `HOMarkov` library

```
matrix([[0.12725632, 0.07761733, 0.07851986, 0.066787  , 0.09566787, 0.05415162, 0.02617329, 0.03790614, 0.10469314, 0.04241877, 0.06588448, 0.03971119, 0.0433213 , 0.03158845, 0.07851986, 0.02978339],
        [0.13055556, 0.08194444, 0.09861111, 0.09722222, 0.11111111, 0.05277778, 0.02222222, 0.05833333, 0.02638889, 0.02638889, 0.03055556, 0.04722222, 0.04027778, 0.04722222, 0.08055556, 0.04861111],
        [0.12808989, 0.07865169, 0.10449438, 0.06404494, 0.07191011, 0.05280899, 0.0247191 , 0.05730337, 0.13033708, 0.04382022, 0.0494382 , 0.04269663, 0.01573034, 0.03483146, 0.05505618, 0.04606742],
        [0.04096045, 0.03672316, 0.06920904, 0.05225989, 0.0819209 , 0.04237288, 0.01694915, 0.05225989, 0.11016949, 0.06920904, 0.16666667, 0.06638418, 0.03531073, 0.04519774, 0.06638418, 0.0480226 ],
        [0.09655938, 0.06659267, 0.07436182, 0.0654828 , 0.0854606 , 0.05216426, 0.03107658, 0.05105438, 0.08324084, 0.06770255, 0.05993341, 0.04106548, 0.04772475, 0.03773585, 0.09655938, 0.04328524],
        [0.12619503, 0.10707457, 0.11281071, 0.10516252, 0.08604207, 0.04397706, 0.02676864, 0.0458891 , 0.04015296, 0.01720841, 0.03632887, 0.0210325 , 0.05927342, 0.05353728, 0.06118547, 0.05736138],
        [0.07662835, 0.06896552, 0.07279693, 0.04980843, 0.07279693, 0.03831418, 0.03065134, 0.05363985, 0.08429119, 0.05363985, 0.07279693, 0.05363985, 0.03831418, 0.0651341 , 0.10727969, 0.06130268],
        [0.04693141, 0.04873646, 0.07942238, 0.06137184, 0.07942238, 0.05415162, 0.0198556 , 0.06498195, 0.0866426 , 0.05595668, 0.12454874, 0.06859206, 0.03429603, 0.05415162, 0.05595668, 0.06498195],
        [0.12807377, 0.07172131, 0.10758197, 0.04713115, 0.08606557, 0.04508197, 0.02254098, 0.04815574, 0.10245902, 0.04508197, 0.07581967, 0.0317623 , 0.03278689, 0.03893443, 0.07377049, 0.04303279],
        [0.102     , 0.072     , 0.106     , 0.058     , 0.118     , 0.044     , 0.044     , 0.052     , 0.028     , 0.03      , 0.032     , 0.026     , 0.082     , 0.044     , 0.112     , 0.05      ],
        [0.17662008, 0.06988564, 0.12452351, 0.07878018, 0.05082592, 0.05209657, 0.01651842, 0.04828463, 0.09656925, 0.0292249 , 0.04193139, 0.03557814, 0.03684879, 0.0317662 , 0.06988564, 0.04066074],
        [0.03353057, 0.04142012, 0.03944773, 0.0295858 , 0.09467456, 0.0433925 , 0.02366864, 0.04733728, 0.08678501, 0.05522682, 0.14792899, 0.0887574 , 0.04733728, 0.0729783 , 0.09467456, 0.05325444],
        [0.07954545, 0.04772727, 0.04772727, 0.05454545, 0.06136364, 0.05681818, 0.03409091, 0.04772727, 0.09772727, 0.07272727, 0.08181818, 0.05      , 0.04090909, 0.06818182, 0.10454545, 0.05454545],
        [0.125     , 0.09475806, 0.08870968, 0.09879032, 0.10483871, 0.04637097, 0.01612903, 0.05846774, 0.03225806, 0.01612903, 0.02419355, 0.02620968, 0.06048387, 0.07459677, 0.08064516, 0.05241935],
        [0.08774038, 0.06490385, 0.046875  , 0.0625    , 0.05528846, 0.03725962, 0.01802885, 0.04927885, 0.16826923, 0.06730769, 0.07692308, 0.07331731, 0.02403846, 0.03966346, 0.07211538, 0.05649038],
        [0.05482042, 0.02646503, 0.03780718, 0.06049149, 0.09829868, 0.05671078, 0.02646503, 0.06994329, 0.09073724, 0.04725898, 0.11153119, 0.05671078, 0.0510397 , 0.06427221, 0.06805293, 0.07939509]])
```

<p align="center">
  <img width="80%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/k%3D2-markov-dengue.png">
</p>

<p align="center">
  <i><b>Figure: </b>Heatmap representation of the transition matrix for k-order markov chain (k = 2)</i>
</p>

#### Compression results using `dmc.c`

`gcc ./library/dmc.c -o dmc`  
`./dmc c ./data/dengue.fa`
```
using 16777216 bytes of predictor memory
compress done: bytes in 10736, bytes out 2876, ratio 0.267884
```
The compression ratio for first order dynamic markov compression is **26.78%**

## Built with
* [Biopython](https://biopython.org) - The Biopython Project is an open-source collection of non-commercial Python tools for computational biology and bioinformatics, created by an international association of developers
* [JupyterLab](https://github.com/jupyterlab/jupyterlab) - Next-generation web-based user interface for Project Jupyter
* [Matplotlib](https://matplotlib.org/) - Plotting library for the Python programming language
* [NumPy](https://www.numpy.org/) - Library for the Python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays
* [SciPy](https://www.scipy.org/) - Free and open-source Python library used for scientific computing and technical computing
* [seaborn](https://seaborn.pydata.org/) - A Python data visualization library based on matplotlib

## Authors

- **Rajat Dipta Biswas** - *Initial work*

See also the list of [contributors](https://github.com/rajatdiptabiswas/dna-compression/graphs/contributors) who participated in this project.

## Acknowledgements

- [NCBI](https://www.ncbi.nlm.nih.gov) | The National Center for Biotechnology Information
- [ROSALIND](http://rosalind.info/glossary/) | An educational resource and web project for learning bioinformatics through problem solving and computer programming
- GitHub | [iz4vve/HOMarkov](https://github.com/iz4vve/HOMarkov)
- Dynamic Markov Compression implementation in C | G. V. Cormack, R. N. S. Horspool | [dmc.c](https://www.jjj.de/crs4/dmc.c)
- YouTube | Google Developers | [Markov Chain Compression (Ep 3, Compressor Head)](https://www.youtube.com/watch?v=05RFEGWNxts)
- YouTube | itechnica | [Dynamic Markov Compression with example](https://www.youtube.com/watch?v=-B_d9RfwI2M)
- Documentations - [Biopython](https://biopython.org/wiki/Documentation), [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/), [Matplotlib](https://matplotlib.org/contents.html), [NumPy](https://docs.scipy.org/doc/numpy-1.13.0/reference/), [SciPy](https://docs.scipy.org/doc/scipy/reference/), [pandas](https://pandas.pydata.org/pandas-docs/stable/), [seaborn](https://seaborn.pydata.org/api.html)

## References

<b><a id="1">[1]</a></b> G. V. Cormack, R. N. S. Horspool, **"Data Compression Using Dynamic Markov Modelling"**, *The Computer Journal*, Volume 30, Issue 6, December 1987, Pages 541–550, https://doi.org/10.1093/comjnl/30.6.541  

**[2]** Wandelt, Sebastian & Bux, Marc & Leser, Ulf, **"Trends in Genome Compression"**, *Current Bioinformatics*, Volume 9, Issue 3, May 2014, https://www.researchgate.net/publication/263474675_Trends_in_Genome_Compression  

**[3]** Whitehead, R. Fletcher. (1994). **"An exploration of dynamic Markov compression"**, https://ir.canterbury.ac.nz/handle/10092/9572  
