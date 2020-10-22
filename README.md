# Shortest Systoles in H(1,1)

This program computes systolic rations of Origamis in the stratum H(1,1). The
algorithm is described in [this preprint](https://arxiv.org/abs/1809.10327v2)
and the code is documented and includes examples.

An example session might look as follows:
    
    sage: load('systoles_h11.pyx')
    Compiling (...)/systoles_h11.pyx...
    sage: max_systole_h11(8)
    ((1,2)(3,4,5,6,7,8)
     (1,3,2,4,6,8)(5,7), (2, 'loop', (1, 0)))
     
## surface_dynamics

The program uses the package `surface_dynamics` by Vincent Delecroix that
available at [GitLab](https://gitlab.com/videlec/surface_dynamics). To install `surface_dynamics` in your installation of Sage, you may issue the command 

    sage -pip install surface_dynamics

in your favourite shell. If your Sage installation is a system-wide
installation and you have no administrator privileges, then you can install the
package `surface_dynamics` for your user only by issuing the command

    sage -pip install surface_dynamics --user

## sage-pylint

The code was checked for adherence to coding standards with `pylint`. You can
install `pylint` by running the command

    sage -pip install pylint
    
in your favourite shell. After a successful installation you should be able to
run the bundled `sage-pylint` script with `sage-pylint systoles_h11.pyx`.


    
