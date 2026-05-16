#!/bin/zsh

# \documentclass{article}
# \usepackage{amsmath,mathrsfs}
# \pagestyle{empty}
# \begin{document}
# \begin{align*}
#   \frac{d \kappa_{xx}(\omega)}{d \omega}
# \end{align*}
# \end{document}

latex $1.tex
dvips -E $1.dvi -o $1.eps
gs -q -dNoOutputFonts -dBATCH -dNOPAUSE -sDEVICE=eps2write -sOutputFile=$1-ol.eps $1.eps
