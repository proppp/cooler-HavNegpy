Introduction
============

``HavNegpy`` is an Python package to analyze dielectric loss and real part of conductivity.
``HavNeg`` is an acronym for Havriliak and Negami function.
 
The package contains three modules to analyze the dielectric data. Each module can be instantiated as::

              >>> import HavNegpy as h
              >>> hn = h.HN()
              >>> hn_deri = h.HN_derivative()
              >>> cond = h.Conductivity()
	   

All modules contain same methods which includes:
``selecting range of data``, ``dumping initial fit parameters``, ``perform least squares fitting``, ``creating an analysis file to save fit results``, and ``save the fit results``.
Besides, other method include ``initial view of fit parameters``

A clear description is provided in the tutorial.
<<<<<<< HEAD
Note: Some HTML images aren't rendered in the tutorials page. All tutorial notebooks are available at `<https://github.com/mkolmang/Tutorials_HavNegpy>`_

=======
Note: Some HTML images aren't rendered in the tutorials page. All tutorial notebooks are available at https://github.com/mkolmang/Tutorials_HavNegpy 
>>>>>>> 0956321fd5ec8a3db3f73a84e2f0420eade0ab49



