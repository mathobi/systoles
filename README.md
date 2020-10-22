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
    
