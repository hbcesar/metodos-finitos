void printTime(GMRES_ParametersPtr par, Solution s){
  //Save the table with each matrix infos in a txt File
  File *log = *fopen( strcat(strcat("testes/", par->name), "_log.txt"), "a" );

  unsigned int kmax = par->kmax[par->krilov_space_index];
  unsigned int fill_in = par->fill_in[par->ilu];
  int reordering = par->reordering;
  int preconditioner = par->preconditioner;


  fprintf(log, "Condicionamento: %d\n", preconditioner);
  fprintf(log, "Reordenamento: %d\n", reordering);
  fprintf(log, "ILU Fill in: %d\n", fill_in);
  fprintf(log, "K: %d\n", kmax);
  fprintf(log, "Iterações: %d\n", s.iterations);
  fprintf(log, "Tempo total: %.4lfs\n\n\n", s.time);

  fprintf(log, "\n\n", );

  fclose(log);
}
