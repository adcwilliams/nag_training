project: 1D heat equation
src_dir: ./src
output_dir: ./doc
summary: Fortran Workshop
author: Andrew Williams
author_description: Fortran workshop attendee
email: a.d.c.williams@bristol.ac.uk
linkedin: 
docmark: !
predocmark: > 
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true
graph: true
coloured_edges: true
search: true
warn: false
license: by-nc
version: 1.1.1

This is the solution for the [Fortran Modernisation Workshop](https://www.nag.co.uk/content/fortran-modernization-workshop).
The workshop exercise solves. 
\begin{equation}
\frac{\partial H}{\partial t} - \kappa\frac{\partial^{2} H}{\partial x^{2}} = f(x)
\end{equation}

@Note
This is a series of Fortran workshops to help computational scientists write
efficient and portable code, and to encourage best practices in software
engineering. Modified by Andrew Williams

@Bug
If you find any bugs, please don't email [me](a.d.c.williams@bristol.ac.uk)
