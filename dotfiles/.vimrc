augroup vimrcEx
  au BufRead * if line("'\"") > 0 && line("'\"") <= line("$") |
  \ exe "normal g`\"" | endif
augroup END 

" An example for a .vimrc file
" Maintener : ken
"
autocmd BufNewFile,BufRead *.py nnoremap <C-e> :!python %
autocmd BufNewFile,BufRead *.rb nnoremap <C-e> :!ruby %
autocmd BufNewFile,BufRead *.pl nnoremap <C-e> :!perl %

set updatecount=0 
set nobackup
set bs=2        " allow backspacing over everything in insert mode
set ai          " always set autoindenting on
"set backup     " keep a backup file
set ruler
set shiftwidth=4
set showmatch
set noautoindent
set updatecount=0 " no swap file"
set nu          " enabling to express number of line
set notitle
set selection=old
set tabstop=4
set viminfo='20,\"50    " read/write a .viminfo file, don't store more
            " than 50 lines of registers
" nobeep
set vb t_vb=
" colour
"#set t_mret t_so=
"#t restore screen
"#set t_ti= t_te=

"syntax enable
colorscheme desert    

call plug#begin('~/.vim/plugged')

Plug 'lervag/vimtex'

call plug#end()

set nocompatible

" Load Vimtex
let &rtp  = '~/.vim/bundle/vimtex,' . &rtp
let &rtp .= ',~/.vim/bundle/vimtex/after'

" Load other plugins, if necessary
" let &rtp = '~/path/to/other/plugin,' . &rtp

filetype plugin indent on
syntax enable

" Vimtex options go here

let g:latex_latexmk_options = '-pdf'
let g:vimtex_compiler_latexmk = {'callback' : 0}

"set spell 
set formatoptions=t
set textwidth=80

set tw=0
