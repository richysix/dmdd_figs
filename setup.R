# setup directories
for ( new_dir in c('data', 'plots', 'output') ) {
  dir_path <- file.path(wd, new_dir)
  if( !dir.exists(dir_path) ){
    dir.create(dir_path)
  }
}
