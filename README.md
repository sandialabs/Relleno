relleno is a SageMath library that facilitates derivations of numerical flux functions for computational solution of conservation laws. 
Its main application has been the derivation of entropy stable flux functions for complex compressible gas flows. 
However, the fundamental library itself is quite general.


# Clone, Install, and Test
To use relleno you'll need SageMath installed on your computer (see http://www.sagemath.org/index.html), and to clone the repo here:

```
git clone https://github.com/sandialabs/relleno.git
```

This will create a directory named `relleno`.
Enter into this directory and install with the following command.

```
cd relleno
sage -pip install --upgrade --no-index -v .
```

To test that the installation was successful, run the following command from the `relleno` directory.
This might take up to thirty seconds or so.

```
sage relleno/tests.sage
```