# setup directories
wd <- getwd()
for ( new_dir in c('plots', 'output') ) {
  dir_path <- file.path(wd, new_dir)
  if( !dir.exists(dir_path) ){
    dir.create(dir_path)
  }
}

file.create('setup.done')
