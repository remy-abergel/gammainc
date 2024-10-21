# Fast and accurate evaluation of a generalized incomplete gamma function

## Purpose

This software focuses on the computation of the (generalized
incomplete Gamma) integral

```math
I_{x,y}^{\mu,p} = \int_{x}^{y} s^{p-1}\,e^{-\mu s} \, ds
```

for $`0 \leq x < y \leq +\infty`$, $`\mu = \pm1`$ and $`p > 0`$. 
Notice that when $`\mu = -1`$, we impose the restriction that $`y \neq
+\infty`$ and $`p`$ is an integer in order to ensure that the integral
is finite and real valued. Details and explanations about the
algorithm are available in the [companion
paper](https://dl.acm.org/doi/10.1145/3365983) of this algorithm (PDF
available [here](https://dl.acm.org/doi/pdf/10.1145/3365983)).
   
**R. Abergel and L. Moisan**. *Algorithm 1006: Fast and Accurate
Evaluation of a Generalized Incomplete Gamma Function*, ACM
Transactions on Mathematical Software, Volume 46, Issue 1, March 2020,
DOI : [10.1145/3365983](https://doi.org/10.1145/3365983).

The *admissible parameter range* for our algorithm is $`0 \leq x,y\leq
10^{15}`$ and $`\varepsilon_{\text{machine}} < p \leq 10^{15}`$ (in
standard double floating point precision) and, in our paper, we show
that the accuracy reached by our algorithm is essentially optimal
considering the limitations imposed by the floating point
arithmetic. This large range of admissible parameters is made possible
by using an appropriate normalization of the generalized incomplete
gamma integral, different from the classical normalization by the
complete Gamma integral. The normalized quantity introduced in our
paper is the *G-function*, which is defined by

```math
\forall x\in\mathbb{R}\cup\{+\infty\}\,,~\forall p > 0\,,\quad G(p,x) = 
\left\{\begin{array}{cl}
\displaystyle{e^{x-p\,\log{(|x|)}} \, \int_{0}^{|x|} s^{p-1} \, e^{-\frac{x}{|x|} s}} \, ds &\text{if } x\leq p\\[7pt]
\displaystyle{e^{x-p\,\log{(x)}} \, \int_{x}^{+\infty} s^{p-1} \, e^{-s}} \, ds&\text{otherwise.}
\end{array}\right.
```

In our paper, we show how the value of $`I_{x,y}^{\mu,p}`$ can be
derived from this normalized quantity which has the avantage to not
being subject to underflow or overflow issues on a large domain.
   
Feel free to use and adapt this code. However, if you use it for your
research or in software code, please be so kind as to cite the paper.

```bibtex
@article{10.1145/3365983,
  author = {Abergel, R\'{e}my and Moisan, Lionel},
  title = {Algorithm 1006: Fast and Accurate Evaluation of a Generalized Incomplete Gamma Function},
  year = {2020},
  issue_date = {March 2020},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {46},
  number = {1},
  issn = {0098-3500},
  url = {https://doi.org/10.1145/3365983},
  doi = {10.1145/3365983},
  journal = {ACM Transactions on Mathematical Software},
  month = {mar},
  articleno = {10},
  numpages = {24}
}
```

## Installation
### Dependancies and content

This software has **very low dependencies**. Only the following
standard libraries are required :

+ `stdio`
+ `stdlib`
+ `math`
+ `string`

This software is written in C language and can be used directly from
the terminal through the provided command line interface. An
additional Matlab interface (MEX) is also provided. The software
content is summarized below.

| source file name             | description |
|------------------------------|-------------|
| [kernel.c](./src/kernel.c)                     | this file contains **public routines** for the evaluation of the G-function and the generalized incomplete gamma function, using either some continued fractions, a recursive integration by parts or a Romberg approximation. |
| [Gfunc.c](./src/Gfunc.c)                      | **command line interface** for the computation of the G-function. The input parameters $`(p,x)`$ must be entered by the user, then, the program displays the value of $`G(p,x)`$ on the standard output with 17 digits of precision. |
| [deltagammainc.c](./src/deltagammainc.c)              | **command line interface** for the computation of the of the generalized incomplete gamma function, $`I_{x,y}^{\mu,p}`$, with a mantissa-exponent representation. The input parameters $`(x,y,\mu,p)`$ must be entered by the user. The computation method for the integral is automatically selected by the algorithm, according to the value of $`(x,y,\mu,p)`$. This program displays into the standard output the following quantities:<br><br>- The computed **mantissa** $`\rho`$, with 17 digits (double precision)<br>- The computed **exponent** $`\sigma`$, with 17 digits (double precision)<br>- The **value of the integral** $`\rho\cdot e^{\sigma}`$, with 17 digits. |
| [deltagammainc_mexinterface.c](./src/deltagammainc_mexinterface.c) | **MEX interface** for the computation of the generalized incomplete gamma function $`I_{x,y}^{\mu,p}`$ from MATLAB. |
| [deltagammainc.m](./src/deltagammainc.m)              | **MATLAB wrapper** calling the MEX interface (after performing consistency checks of the inputs) and performing the (string) base 10 scientific notation formatting of the computed integral. |
| [INSTALL](./src/INSTALL)                      | executable file for C-module **compilation** (bash) |
| [kernel.h](./src/kernel.h)                     | **header** file |
				   
### Compilation of the C modules

See below how the program can be compiled using gcc (the compilation
commands can be adapted to others C-compiler).

You can compile the program using the [INSTALL](src/INSTALL)
executable. To that aim, open a terminal, place yourself into the
`src` directory of this program, and simply run the command

```bash
./INSTALL
```

When the compilation is successful, the following message is displayed
in the standard output:

```
*************************************************************
*           COMPILATION OF DELTAGAMMAINC MODULES            *
*************************************************************
	
+ compilation of module 'Gfunc': success
+ compilation of module 'deltagammainc': success

```

Notice that you can also manually compile the program using the
following commands

```bash
gcc -O3 kernel.c Gfunc.c -lm -o Gfunc
gcc -O3 kernel.c deltagammainc.c -lm -o deltagammainc
```

When the compilation is done, the executable files `deltagammainc` and
`Gfunc` are created into the `src` directory of this program. Those
executable files are command line interfaces that allows to evaluate
the G-function $`G(p,x)`$ and the generalized incomplete Gamma function
$`I_{x,y}^{\mu,p}`$. 

**Display program documentation**

Executed from the `src` directory, the commands 

```bash
./Gfunc --help
./deltagammainc --help
```

will display the documentation of the `Gfunc` and `deltagammainc`
programs, that is 

```
Usage: Gfunc [--help] x p

   --help : display help
   x      : (double) a real number, possibly infinite (-inf < x <= inf)
   p      : (double) a positive real number (p > 0)

Description: compute G(p,x) such as

  if x <= p: G(p,x) = exp(x-p*log(|x|)) * integral over [0,|x|] of s^{p-1} * exp(-sign(x)*s) ds
  otherwise: G(p,x) =  exp(x-p*log(x))  * integral over [x,inf] of s^{p-1} * exp(-s) ds

```

and 

```
Usage: deltagammainc [--help] [--verbose] x y mu p

   --help    : display help
   --verbose : verbose mode (display how the integral was computed)
   x         : (double) a nonnegative number, possibly infinite (inf)
   y         : (double) a number, possibly infinite (inf), greater than x
   mu        : (double) a nonzero real number
   p         : (double) a positive real number (must be integer if mu<0)

Description:

   Compute (rho,sigma) such as I = rho * exp(sigma), where
   I = integral over [x,y] of s^{p-1} * exp(-mu*s) ds.

```

**Commented basic example**
	
Let us compute the integral $`I_{x,y}^{\mu,p}`$ for $`x=0.1`$,
$`y=200`$, $`p=100`$ and $`\mu=1`$. From the `src` directory, run
the following command:

```bash
./deltagammainc --verbose 0.1 200 1 100
```

This should result in the following displayed information (without the line
numerotation, that were added here to improve the quality of the
explanations)

```
1.  Computing I = integral over [x,y] of s^{p-1} * exp(-mu*s) ds,
2.  where (hexadecimal): x=0x1.999999999999ap-4, y=0x1.9p+7, mu=0x1p+0, p=0x1.9p+6,
3.  (in decimal approx): x=0.1000000000000000055511151231257827, y=200, mu=1, p=100.
4.
5.  Details: I will be computed as a difference of integrals I=A-B.
6.
7.  Output: I = 9.33262154439491644e+155
8.
9.  Representation with double numbers: I = rho*exp(sigma), where
10.  rho   = 9.99999999999998113e-01
11.  sigma = 3.59134205369575454e+02
12.
```

Let us comment this result line by line:

+ **line 1**: we recall the definition of $`I_{x,y}^{\mu,p}`$

+ **line 2**: we display the values of the inputs $`(x,y,mu,p)`$ in
  hexadecimal floating-point format. The interest of using this format
  is that all hexadecimal floating-point constants have exact
  representations in binary floating-point, unlike decimal
  floating-point constants, which in general do not (for further
  details on the syntax of hexadecimal floating-point constants, see
  pages 57-58 of the C99 specification).

+ **line 3**: decimal approximation of the inputs $`(x,y,mu,p)`$. When
  an input value is not exactly representable in double floating-point
  precision (such as $`x=0.1`$), it is replaced by the nearest double
  floating-point number (here, the actual value of $`x`$ will be
  approximatively equal to `0.1000000000000000055511151231257827`, see
  line 2. for the exact representation of x in hexadecimal
  floating-point format).

+ **lines 4-6**: because the option --verbose was used, we can read
  here wether $`I_{x,y}^{\mu,p}`$ was numerically computed using a
  difference or using Romberg's numerical integration method (we refer
  to the companion paper for more details about the effective
  numerical procedure).

+ **line 7**: we display the computed value of $`I_{x,y}^{\mu,p}`$, in
  scientific notation.

+ **lines 8-12**: $`I_{x,y}^{\mu,p}`$ was first computed with a
  mantissa-exponent representation of the type $`I_{x,y}^{\mu,p} =
  \text{rho}\cdot \exp(\text{sigma})`$, where rho and sigma where
  evaluated in double precision (their computed values are displayed
  in scientific notation with 17 digits of precision).

### Compilation of the MEX files (for Matlab interfacing only)

To compile the MEX interface, open a Matlab console and place yourself
into the `src` directory of this software. Then, execute the
following command : 

```matlab
mex -R2018a -silent -lm CFLAGS="\$CFLAGS -std=c99" kernel.c deltagammainc_mexinterface.c -output deltagammainc_mexinterface
```

A warning message related to unsupported gcc compiler may be raised
and can be ignored. When the compilation is successful, the
`deltagammainc_mexinterface.mexa64` file is created in the `src`
directory of the software. The `deltagammainc_mexinterface` command
can then be called from the Matlab console. For instance, the integral
$`I_{x,y}^{\mu,p}`$ for $`(x=0.1,y=200,\mu = 1, p = 100)`$ can be
evaluated using the command

```matlab
x = 0.1; y = 200; mu = 1; p = 100; 
[rho,sig] = deltagammainc_mexinterface(x,y,mu,p); 
fprintf("rho = %.17e\n",rho); 
fprintf("sigma = %.17e\n",sig); 
```

which should display the following message into the Matlab console:

```
rho = 9.99999999999998113e-01
sigma = 3.59134205369575511e+02
```

where the values of $`\rho`$ and $`\sigma`$ correspond to the
mantissa-exponent representation of the integral and satisfy
$`I_{x,y}^{\mu,p} = \rho \cdot e^{\sigma}`$.

Instead of using `deltagammainc_mexinterface`, we recommend the use of
the provided Matlab wrapper `deltagammainc.m` which performs
consistency checks of the inputs and format the mantissa-exponent
representation of the integral into (string) base 10 scientific
notation. Indeed, the command

```matlab
x = 0.1; y = 200; mu = 1; p = 100; 
[rho,sig,strval] = deltagammainc(x,y,mu,p);
disp(strval); 

```

displays the computed value of $`I_{x,y}^{\mu,p}`$ (for $`x=0.1, y=200, \mu = 1, p=100`$):

```
9.33262154439544803e+155
```

Evaluation of $I_{x,y}^{\mu,p}$ for multiple values of $(x,y,\mu,p)$
is also possible. When $`(x,y,\mu,p)`$ are arrays containing $N$ elements each, the `deltagammainc` module returns three arrays `rho`,`sigma` and `strval` (with same shape as `x`) such as `rho(k)`, `sigma(k)` and `strval(k)` correspond to the computed mantissa, exponent and formatted value of the integral $`I_{x(k),y(k)}^{\mu(k),p(k)}`$ for $1 \leq k \leq N$. For instance, running the commands

```matlab
x = [0.1,10]; y = [100,200]; mu = [-1,1]; p = [50,300]; 
[rho,sig,strval] = deltagammainc(x,y,mu,p);
disp(strval);
```

yields 

```
"1.80009633885347035e+141"    "2.76621762074220223e+601"
```

which correspond to the computed values of the integral
$`I_{x,y}^{\mu,p}`$ for $`(x,y,\mu,p) = (0.1,100,-1,50)`$ and
$`(x,y,\mu,p) = (10,200,1,300)`$.

## Licensing

  DELTAGAMMAINC Fast and Accurate Evaluation of a Generalized
  Incomplete Gamma Function. Copyright (C) 2016 Remy Abergel
  (remy.abergel AT u-paris.fr), Lionel Moisan (Lionel.Moisan AT
  u-paris.fr).

  DELTAGAMMAINC is a software dedicated to the computation of a
  generalized incomplete gammafunction. See the Companion paper for a
  complete description of the algorithm:

  R. Abergel and L. Moisan. ``Algorithm 1006: Fast and Accurate
  Evaluation of a Generalized Incomplete Gamma Function'', ACM
  Transactions on Mathematical Software, Volume 46, Issue 1, March
  2020, DOI : <https://doi.org/10.1145/3365983>.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

## Contact

  This work was done at [Université Paris Cité](https://u-paris.fr/),
  [Laboratoire MAP5 (CNRS UMR
  8145)](https://map5.mi.parisdescartes.fr/), 45 rue des Saints-Pères
  75270 Paris Cedex 06, FRANCE.

  If you have any comments, questions, or suggestions regarding this
  code, or if you find a bug, don't hesitate to open a
  [discussion](https://github.com/remy-abergel/gammainc/discussions)
  or a [bug issue](https://github.com/remy-abergel/gammainc/issues).

  Also, if you find this code useful, we would be delighted to now
  about what application you are using it for.

  Rémy Abergel & Lionel Moisan.
