`Relleno` is a SageMath library that facilitates derivations of numerical flux functions for computational solution of conservation laws. 
Its main application has been the derivation of entropy stable flux functions for compressible gas flows. 
However, the fundamental library itself is quite general.


# Clone, Install, and Test
To use Relleno you'll need SageMath installed on your computer (see http://www.sagemath.org/index.html).
After installing SageMath, clone the `Relleno` repository:

```
git clone https://github.com/sandialabs/Relleno.git
```

This will create a directory named `Relleno`.
Enter into this directory and install with the following command.

```
cd Relleno
sage -pip install --upgrade --no-index -v .
```

To test that the installation was successful, run the following command from the `Relleno` directory.
This might take up to thirty seconds or so.

```
sage Relleno/tests.sage
```

# Theory Manual
The `docs/Relleno.pdf` file details the theory underlying `Relleno`. It probably won't be too helpful for learning how to use the code, though.

# Demonstrations
To see `Relleno` in use, check out the `demo` folder.
There are several `.ipynb` files for use with the Jupyter notebook interface to SageMath (`sage -n jupyter`).
