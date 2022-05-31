from numpy.distutils.core import Extension

ext = Extension(name = 'fortran',
                 sources = ['Legendre_evals.f90', 'Integral_evals.f90'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'fortran_for_pygeodyn',
          description       = "F2PY Users Guide examples",
          author            = "Pearu Peterson",
          author_email      = "pearu@cens.ioc.ee",
          ext_modules = [ext]
          )
# End of setup_example.py
