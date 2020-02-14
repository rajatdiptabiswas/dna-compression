# DNA Compression Research

Analyzing compression algorithms for genomic sequencing data

The  Jupyter Notebooks in the repository contain further information

- **`dna-notes.ipynb`** - Notes made from different research papers and websites.
- **`dna-sample-analysis.ipynb`** - Analyzes the protein sequence in `data/sample.fa` and gives information about the frequency of different nucleobases in the sample FASTA sequence.
- **`dmc-compress-test.ipynb`** - Compresses `data/human_data.txt` to `data/human_data_compressed` and then decompresses it back to `data/human_data_uncompressed.txt` using `library/dmc.c`. Gives implementation details on how to compress and decompress files using the C implementation used in the paper[<sup>[1]</sup>](#1). Calculates the compression ratio.
- **`dynamic-markov-compression.ipynb`** - Implements ***dynamic Markov compression*** algorithm proposed in the paper by G. V. Cormack and R. N. Horspool [<sup>[1]</sup>](#1) to find the compression ratio on the genomic data of a dengue virus `data/dengue.fa`. Visualizations made using the `HOMarkov` library.
- **`markov-arithmetic-implementation.ipynb`** - Custom personal implementation of a k-th order Markov chain to compress a string using arithmetic coding.

**Dynamic Markov Compression** (DMC) uses predictive arithmetic coding similar to the prediction by partial matching (PPM), except that the input is predicted one bit at a time (rather than one byte at a time). DMC has a good compression ratio and moderate speed, similar to PPM, but requires somewhat more memory and is not widely implemented.

Implementing ***Markov modeling*** and ***Guazzo's arithmetic coding*** together provides a powerful data compression method. Its advantage is that it is adaptive, i.e. messages can be encoded and decoded with a single pass through the data.

All compression algorithms rely on a priori assumption of the data. The new direction taken here is that an algorithmic attempt is used to discover a Markov chain model that describes the data. If it can be created using the first part of the data then it can be used to predict forthcoming characters.

The decoder will generate the same statistical properties that the encoder is using for the data as the decoder is decoding it.

## Theory

Data compression can be divided into 2 parts
- **Modeling** - Statistical Modeling, Markov Modeling (aka Finite Context Modeling)
- **Encoding** - Huffman Coding, Arithmetic Coding

**Modeling** is done to calculate probabilities for the encoding process.  
**Encoding** does the actual compression. Fewer bits are used to determine the same data.

**Markov chains** are mathematical systems that experience state transitions based on probabilistic rules.
They can simply be defined as a Finite State Machine with probabilities on the edges of the nodes of the graph.

A k-order Markov chain requires <b>4<sup>m</sup> states</b>. Hence the chains get memory intensive very quickly.

First-order Markov chains look like the following. They have 4<sup>1</sup> = 4 states.

<p align="center">
  <img width="50%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/first-order-chain.png">
</p>
<p align="center">
  <i>First-order Markov Chain example</i>
</p>

Second order markov chains require 4<sup>2</sup> = 16 states.

<p align="center">
  <img width="60%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/second-order-chain.png">
</p>
<p align="center">
  <i>Second-order Markov Chain example</i>
</p>

This can be represented in a 16 x 4 table because all the states will need 4 transitions.

<table align="center"><tr><th></th><th>P(A | row)</th><th>P(T | row)</th><th>P(C | row)</th><th>P(G | row)</th></tr><tr><td>AA</td><td> </td><td> </td><td> </td><td> </td></tr><tr><td>AT</td><td> </td><td> </td><td> </td><td> </td></tr><tr><td>...</td><td> </td><td> </td><td> </td><td> </td></tr><tr><td>GG</td><td> </td><td> </td><td> </td><td> </td></tr></table>

The figure below shows how finite-context models are implemented. The rows of the table represent probability models at a given instant t. In this example, the particular model that is chosen for encoding a symbol depends on the last five encoded symbols (order-5 context).

<p align="center">
  <img width="50%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/probability-table.png">
</p>
<p align="center">
  <i>Simple example illustrating how finite-context models are implemented</i>
</p>

### Higher-order context models
Contexts of 8th order tend to work best in most cases, better than everything except infinite order.

While order-8 models achieve the best compression of the data itself, it's often been found that anything beyond order-5 is not worth it.

### Cloning
In actual models, all the states are not made in the implementation. Cloning is used.
States are cloned to ‘remember’ previous states. If enough transitions occur (figured out using threshold values) through the new cloned state, the cloned state is kept in the final model.

<p align="center">
  <img width="60%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/cloning.png">
</p>
<p align="center">
  <i><b>(a)</b> Part of a Markov model; <b>(b)</b> the Markov model after "cloning"</i>
</p>

### Entropy
For lossless compression there is a limit to how much one can compress the data. This is calculated using entropy.

<p align="center">
  <img width="50%" src="https://render.githubusercontent.com/render/math?math=Entropy%5C%20H(X)%20%3D%20-%5Csum%20p(X)%5Clog%20p(X)">
</p>

For equally probable 4 items the entropy is 2 bits. It cannot be lower than that. 00, 01, 10 and 11.

However, if the probabilities are skewed for any one of the items, one can have fewer bytes than the other. This leads to savings.

If the probability is high, the formula will give a lower value. 
It is multiplied with a probability to give it a weight.
We find the weighted average of all the characters in the string.


### Modes of compression

#### Static
- Calculate probability
- Generate variable length codes
- Encode

A static model is a fixed model that is known by both the compressor and the decompressor and does not depend on the data that is being compressed. For example, the frequencies of symbols in the English language.

#### Semi-static
A semi-adaptive or semi-static model is a fixed model that is constructed from the data to be compressed. For example, the symbol frequencies computed from the text to be compressed can be used as the model. The model has to be included as a part of the compressed data.

#### Adaptive/Dynamic
- Read symbol
- Generate current variable length codes
- Update variable length codes

An adaptive model changes during the compression. At a given point in compression, the model is a function of the previously compressed part of the data. Since that part of the data is available to the decompressor at the corresponding point in decompression, there is no need to store the model.

The only two things that matter are making sure that 
1) the model attempts to accurately predict the probability when a character will appear, and 
2) the encoder and decoder have identical models at all times.

<p align="center">
  <img width="95%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/frequency-calculation.png">
</p>
<p align="center">
  <i>Symbol frequency calculation example</i>
</p>

### Calculate probabilities
First a frequency table needs to be calculated. For each context, the probabilities are calculated using the formula. (n<sub>A</sub> is the frequency of A in the context)

<b>P(A | context) = (n<sub>A</sub> + c) / (n<sub>A</sub> + n<sub>T</sub> + n<sub>C</sub> + n<sub>G</sub> + 4c)</b>

0 probabilities cause problems further down the line. For this fact, a value c is used. The parameter c controls how much probability is assigned to unseen (but possible) events, and plays a key role in the case of high-order models.

### Encoding
Once the probabilities are calculated, they are used to encode the string to compress.

In Huffman coding, the smaller bits are assigned to higher probable symbols.

In arithmetic coding (the most widely used compression these days), a single decimal number in the range of 0 to 1 is assigned to the whole string. The longer the string the more number of digits after the decimal.

The decimal value is then required to be efficiently stored in computer memory. This is a different research topic in itself. 

#### Arithmetic Coding

Arithmetic coding is a data compression technique that encodes data by creating a code string which represents a fractional value on the number line between 0 and 1. The coding algorithm is symbolwise recursive; i.e. it operates upon and encodes (decodes) one data symbol per iteration or recursion. 

On each recursion, the algorithm successively partitions an interval of the number line between 0 and 1, and retains one of the partitions as the new interval. Thus, the algorithm successively deals with smaller intervals, and the code string, viewed as a magnitude, lies in each of the nested intervals. 

The data string is recovered by using magnitude comparisons on the code string to recreate how the encoder must have successively partitioned and retained each nested subinterval.

<p align="center">
  <img width="95%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/arithmetic-coding.png">
</p>
<p align="center">
  <i>Representation of the Arithmetic Coding Process with the interval scaled up at each stage</i>
</p>

## Getting started

- *(Optional)* Install [Anaconda](https://www.anaconda.com/distribution/#download-section)
- Install [JupyterLab](https://jupyter.org/install)
- Clone the repository
- Create a virtual environment using `requirements/requirements.txt` or `requirements/conda-requirements.txt`
- Open the notebooks using JupyterLab

## Results 

### Contents of the sample file 

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

### Markov chain model created using the `HOMarkov` library

#### k = 1

```
matrix([[0.32340922, 0.21015762, 0.25977817, 0.20665499],
        [0.40223214, 0.23348214, 0.11651786, 0.24776786],
        [0.35234657, 0.18050542, 0.28411552, 0.18303249],
        [0.19147084, 0.21627502, 0.36205396, 0.23020017]])
```
<p align="center">
  <img width="80%" src="https://github.com/rajatdiptabiswas/dna-compression/blob/master/images/k%3D1-markov-dengue.png">
</p>

<p align="center">
  <i><b>Figure 1: </b>Heatmap representation of the transition matrix for k-th order markov chain (k = 1)</i>
</p>


#### k = 2

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
  <i><b>Figure 2: </b>Heatmap representation of the transition matrix for k-th order markov chain (k = 2)</i>
</p>

### Compression results using `dmc.c`

`gcc ./library/dmc.c -o dmc`  
`./dmc c ./data/dengue.fa`
```
using 16777216 bytes of predictor memory
compress done: bytes in 10736, bytes out 2876, ratio 0.267884
```
The compression ratio for dynamic markov compression is **26.78%**

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

**[4]** Langdon, Glen G. (1984). **"An introduction to arithmetic coding"**, *IBM Journal of Research and Development*, https://ieeexplore.ieee.org/abstract/document/5390377

**[5]** P. Krishnamachari, **"The Joy of Finite Context Modeling"**, http://chiranjivi.tripod.com/Finite.html

**[6]** D. Phong, **"Finite Context Modelling"**, http://www.hugi.scene.org/online/coding/hugi%2019%20-%20cofinite.htm

**[7]** Ian & H, Ian & Neal, & M, Radford & Cleary, & G, John. (1987). **"Arithmetic Coding for Data Compression"**, *Communications of the ACM*, https://web.stanford.edu/class/ee398a/handouts/papers/WittenACM87ArithmCoding.pdf

**[8]** M. Mahoney, **"Data Compression Explained"**, http://mattmahoney.net/dc/dce.html

**[9]** **"Text Compression"**, https://www.cs.helsinki.fi/u/tpkarkka/opetus/12k/dct/lecture05.pdf

**[10]** Pinho, Armando & Neves, António & Martins, Daniel & Bastos, Carlos & Ferreira, Paulo. (2010). **"Finite-Context Models for DNA Coding"**, https://pdfs.semanticscholar.org/905b/dfcdf93bdd511f7dcaded65c7aed07da7cef.pdf

**[11]** M. Nelson. (2014). **"Data Compression With Arithmetic Coding"**, https://marknelson.us/posts/2014/10/19/data-compression-with-arithmetic-coding.html

**[12]** M. Nelson. (1991). **"Arithmetic Coding + Statistical Modeling = Data Compression"**, https://marknelson.us/posts/1991/02/01/arithmetic-coding-statistical-modeling-data-compression.html

**[13]** Howard, Paul & Vitter, Jeffrey. (1994). **"Practical Implementations of Arithmetic Coding"**, https://www.cc.gatech.edu/~jarek/courses/7491/Arithmetic2.pdf

**[14]** GeeksforGeeks, **"Floating Point Representation"**, https://www.geeksforgeeks.org/floating-point-representation-basics/

**[15]** Cprogramming.com, **"Floating point number representation"**, https://www.cprogramming.com/tutorial/floating_point/understanding_floating_point_representation.html
