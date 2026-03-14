# Setting luatex as pdf-generating mode

$pdf_mode = 4;
$lualatex = 'lualatex -synctex=1 -interaction=nonstopmode -shell-escape %O %S';

# Work in directory of main tex file

$do_cd = 1;

# Directories for output and aux-files (relative to main tex file)

$out_dir = './';
$aux_dir = 'aux_folder';
ensure_path('TEXINPUTS', './');

# max repeat (default is 5)

$max_repeat=5;

# clean up generated files

$clean_ext = 'aux bbl blg idx lof log lot out toc synctex.gz';