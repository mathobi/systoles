# Shortest Systoles in H(1,1)

This program computes systolic rations of Origamis in the stratum H(1,1). The
algorithm is described in [this preprint](https://arxiv.org/abs/1809.10327v2)
and the code is documented and includes examples. The code was tested with Sage
version 8.9 from 2019-09-29.

A session might look as 
    
    sage: load('systoles_h11.pyx')
    Compiling [...]/systoles_h11.pyx...
    sage: max_systole_h11(8)
    ((1,2)(3,4,5,6,7,8)
     (1,3,2,4,6,8)(5,7), (2, 'loop', (1, 0)))
     
or as

    sage: load('systoles_h11.pyx')
    Compiling [...]/systoles_h11.pyx...
    sage: max_systole_h11(16)
    ((1,2,3)(4,5,6,7,8,9,10,11,12,13,14,15,16)
     (1,4,11,8,5,12)(2,15,9,6,13)(3,16,10,7,14),
     (sqrt(2) + 1, 'cycle', (1, (0, 1)), (sqrt(2), (1, 1))))


     
## surface_dynamics

The program uses the package `surface_dynamics` by Vincent Delecroix whose
source code is available at
[GitLab](https://gitlab.com/videlec/surface_dynamics). You may also have a look
at the [PyPI page](https://pypi.org/project/surface-dynamics) of the package
that contains a lot more of information.

To add the package
`surface_dynamics` to your installation of Sage, you may issue the command 
    sage -pip install surface_dynamics

in your favourite shell. If your Sage installation is a system-wide
installation and you have no administrator privileges, then you can install the
package `surface_dynamics` for your user only by issuing the command

    sage -pip install surface_dynamics --user
    

## sage-pylint

The code was checked for adherence to Python coding standards with `pylint`. You can
install `pylint` by running the command

    sage -pip install pylint
    
in your favourite shell. After a successful installation you should be able to
run the bundled `sage-pylint` script by issuing the command 
`sage-pylint systoles_h11.pyx`. The script `sage-pylint` and hence this step
are *not* necessary to use the program, though.


    
