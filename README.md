# libhxnd

An open source framework for n-dimensional hypercomplex calculations for NMR.

## Introduction

The libhxnd package implements a portable, open-source programming interface
for manipulating a class of hypercomplex numbers found in the field of
[Nuclear Magnetic Resonance] (
http://en.wikipedia.org/wiki/Nuclear_magnetic_resonance) spectroscopy.
More importantly, libhxnd supports the manipulation of these hypercomplex
numbers inside multidimensional arrays of any size. Using libhxnd grammars,
_any_ NMR dataset of _any_ dimensionality may be represented by the same
data structure (`hx_array`).

A few core terms are used when discussing the use of hypercomplex
multidimensional arrays in libhxnd, as follows...

### Algebraic Dimensionality

The algebraic dimensionality of an array (or scalar) refers to the number of
basis _phases_ used to represent the hypercomplex number. For example, a real
number is _zero_-dimensional, because it contains no complex basis phases. A
complex number is _one_-dimensional and contains a single basis phase:

> x = a + b **u1**

On the other hand, a _two_-dimensional hypercomplex number contains two bases,
**u1** and **u2**. A generic two-dimensional hypercomplex number must then
be represented by _four_ real coefficients, like so:

> x = a + b **u1** + c **u2** + d **u1** **u2**

A generic three-dimensional hypercomplex number contains three bases, **u1**,
**u2** and **u3**, and will require _eight_ real coefficients:

> x = a + b **u1** + c **u2** + d **u1** **u2** + e **u3** + f **u1** **u3** +
g **u2** **u3** + h **u1** **u2** **u3**

More generally, the number of coefficients required to express a given
hypercomplex number in _d_ dimensions is two to the power _d_.

### Topological Dimensionality

The topological dimensionality of an array refers to the number of array
dimensions, similar to the concept of tensor order. For example, scalars and
vectors are one-dimensional, matrices are two-dimensional, rectangular cuboids
are three-dimensional, _etc_. A _k_-dimensional array will require a set of
_k_ indices to identify a unique scalar element.

### Topological Size

The topological size of an array refers to the set of sizes of the dimensions
of the array. An array with _k_ topological dimensions will have _k_ sizes,
one for each dimension, _e.g._:

> _k_ = 1: (3)

> _k_ = 2: (4, 5)

> _k_ = 4: (2, 5, 11, 7)

Of course, the size of any dimension may just as easily be **one**, making
that dimension effectively meaningless in the array. For example, an array
of size _(1,1,2,1,8,1)_ is effectively the same as an array of size _(2,8)_.

### Configuration

The configuration of an array refers to the combination of the algebraic
dimensionality, the topological dimensionality, and the topological size
of the array. Only arrays having identical configurations may be added or
multiplied element-wise without modification. Arrays having identical
lengths (total number of elements) may be reshaped for element-wise
arithmetic.

## Representation of NMR data

Because the `hx_array` data structure natively describes NMR data of any
dimensionality, its use in representing such data files requires relatively
few extra pieces of information. Thus, the `datum` structure in libhxnd is
a fairly thin wrapper around the `hx_array` structure, holding per-dimension
metadata such as spectral width, offset and carrier.

## More information

The inspiration for libhxnd arose from the following publication:

> Schuyler A. D., Maciejewski M. W., Stern A. S., Hoch J. C., _Formalism
> for hypercomplex multidimensional NMR employing partial-component
> subsampling_, Journal of Magnetic Resonance, February 2013, 227: 20--24.

## Licensing

The libhxnd project is released under the [GNU GPL 2.0] (
http://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

Enjoy!

*~ Brad.*

