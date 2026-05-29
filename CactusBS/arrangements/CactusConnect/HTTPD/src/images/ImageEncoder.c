/***********************************************************
ImageEncoder: Standalone program to convert any image into
        a static character array (each byte represented in
        hexadecimal) in a header file.  Use this for any
        images that you want to embed in the code so
        that it can run standalone (ie. when you are 
        disconnected from the network).  Currently used
        to spew out the wwwcactuscodeorg.jpg image used
        in thorn_HTTPD through the WebImage.cc extension
        to the HTTPD.
Compilation:
        $(CC) -o ImageEncoder ImageEncoder.c
Running:
        ImageEncoder <inputimage> <outputheader>
        This will take the input binary file and create a
        header file (something.h) which has the image 
        encoded as a static character array with the
        same name as the input image file (minus the extension.
**************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc,char *argv[]){
  char *s;
  char *varname=0;
  int c,counter=0;
  FILE *in,*out;
  if(argc<2) {
    fprintf(stderr,"Usage: %s <inputimage> <outputheader>\n",argv[0]);
    fprintf(stderr,"\tor %s <inputimage>   output to stdout\n",argv[0]);
    exit(0);
  }
  in = fopen(argv[1],"r");
  if(argc<3)
    out=stdout;
  else
    out = fopen(argv[2],"w");
  s = strrchr(argv[1],'/');
  if(s){
    varname = (char*)malloc(strlen(s)+1);
    strcpy(varname,s); /* strip */
  }
  else{
    varname = (char*)malloc(strlen(argv[1])+1);
    strcpy(varname,argv[1]);
  }
  /* printf("VarName=[%s]\n",varname); */
  s = strrchr(varname,'.'); 
  if(s){ /* strip .<ext> */
    *s='\0';
  }
  
  /* Now actually start writing stuff */
  fprintf(out,"#ifndef __%s_H_\n#define __%s_H_\n",
          varname,varname);
  fprintf(out,"\n\n");
  fprintf(out,"static unsigned char %s[]={\n",varname);
  /* inefficient, but it works */
  c=fgetc(in);
  while(c!=EOF){
    int tmp=c;
    c=fgetc(in);
    if(c!=EOF)
      fprintf(out,"0x%02X,",tmp);
    else
      fprintf(out,"0x%02X};",tmp);
    counter++;
    if(!(counter%12)) 
      fprintf(out,"\n\t");
  }
  fprintf(out,"\n\n#endif /* __%s_H_ */\n",varname);
  fclose(in);
  if(argc>2) 
    fclose(out);
  if(varname) free(varname);
}
