/*   Dynamic Markov Compression (DMC)    Version 0.0.0
 
 
     Copyright 1993, 1987
 
     Gordon V. Cormack
     University of Waterloo
     cormack@uwaterloo.ca
 
 
     All rights reserved.
 
     This code and the algorithms herein are the property of Gordon V. Cormack.
 
     Neither the code nor any algorithm herein may be included in any software,
     device, or process which is sold, exchanged for profit, or for which a
     licence or royalty fee is charged.
 
     Permission is granted to use this code for educational, research, or
     commercial purposes, provided this notice is included, and provided this
     code is not used as described in the above paragraph.
 
*/

/*    This program implements DMC as described in

      "Data Compression using Dynamic Markov Modelling",
      by Gordon Cormack and Nigel Horspool
      in Computer Journal 30:6 (December 1987)

      It uses floating point so it isn't fast.  Converting to fixed point
      isn't too difficult.

      comp() and decomp() implement Guazzo's version of arithmetic coding.

      pinit(), predict(), and pupdate() are the DMC predictor.

      pflush() reinitializes the DMC table and reclaims space

      preset() starts the DMC predictor at its start state, but doesn't
               reinitialize the tables.  This is used for packetized
               communications, but not here.

*/

/*
      NOTE: Edited a few lines before compilation

      Rajat D. Biswas
      rdbiswas97@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>


float predict();
int pinit();
int pupdate();

int memsize = 0x1000000;

main(argc,argv)
   int argc;
   char *argv[];
{
   if (argc == 3 && isdigit(*argv[2])) sscanf(argv[2],"%d",&memsize);
   if (argc >= 2 && *argv[1] == 'c') comp();
   else if (argc >= 2 && *argv[1] == 'd') decomp();
   else {
      fprintf(stderr,"usage:  dmc [cd] memsize <infile >outfile\n");
      exit(1);
   }
   return 0;
}

comp(){
   int max = 0x1000000,
       min = 0,
       mid,
       c,i,
       inbytes = 0, 
       outbytes =3,
       pout = 3,
       bit;
   
   pinit(memsize);
   
   for(;;){
      c = getchar();
      if (c == EOF) {
         min = max-1;
         fprintf(stderr,"compress done: bytes in %d, bytes out %d, ratio %f\n",
                         inbytes,outbytes,(float)outbytes/inbytes);
         break;
      }
      for (i=0;i<8;i++){
         bit = (c << i) & 0x80;
         mid = min + (max-min-1) * predict();
         pupdate(bit != 0);
         if (mid == min) mid++;
         if (mid == (max-1)) mid--;
   
         if (bit) { 
            min = mid;
         } else {
            max = mid;
         }
         while ((max-min) < 256) {
            if(bit)max--;
            putchar(min >> 16);
            outbytes++;
            min = (min << 8) & 0xffff00;
            max = ((max << 8) & 0xffff00 ) ;
            if (min >= max) max = 0x1000000;
         }
      }
      if(!(++inbytes & 0xff)){
         if(!(inbytes & 0xffff)){
               fprintf(stderr,
                       "compressing... bytes in %d, bytes out %d, ratio %f\r",
                       inbytes,outbytes,(float)outbytes/inbytes);
         }
         if (outbytes - pout > 256) { /* compression failing */
            pflush();
         }
         pout = outbytes;
      }
   }
   putchar(min>>16);
   putchar((min>>8) & 0xff);
   putchar(min & 0x00ff);
}


decomp(){
   int max = 0x1000000,
       min = 0,
       mid,
       val,
       i,
       inbytes=3,
       pin=3,
       outbytes=0,
       bit,
       c;
   
   pinit(memsize);
   
   val = getchar()<<16;
   val += getchar()<<8;
   val += getchar();
   while(1) {
      c = 0;
      if (val == (max-1)) {
         fprintf(stderr,"expand: input %d output %d\n",inbytes,outbytes);
         break;
      }
      for (i=0;i<8;i++){
         mid = min + (max-min-1)*predict();
         if (mid == min) mid++;
         if (mid == (max-1)) mid--;
         if (val >= mid) {
            bit = 1;
            min = mid;
         } else {
            bit = 0;
            max = mid;
         }
         pupdate(bit != 0);
         c = c + c + bit;
         while ((max-min) < 256) {
            if(bit)max--;
            inbytes++;
            val = (val << 8) & 0xffff00 | (getchar()& 0xff);
            min = (min << 8) & 0xffff00;
            max = ((max << 8) & 0xffff00 ) ;
            if (min >= max) max = 0x1000000;
         }
      }
      putchar(c);
      if(!(++outbytes & 0xff)){
         if (inbytes - pin > 256) { /* compression was failing */
            pflush();
         }
         pin = inbytes;
      }
   }
}

typedef struct nnn {
           float count[2];
           struct nnn    *next[2];
} node;

static int threshold = 2, bigthresh = 2; 

static node *p, *new, nodes[256][256];

static node *nodebuf;
static node *nodemaxp;
static node *nodesptr;

pinit(memsize)
   int memsize;
{
   fprintf(stderr,"using %d bytes of predictor memory\n",memsize);
   nodebuf = (node *) malloc (memsize);
   if (nodebuf == (node *) NULL) {
      fprintf(stderr,"memory alloc failed; try smaller predictor memory\n");
      exit(1);
   }
   nodemaxp = nodebuf + (memsize/sizeof(node)) - 20;
   pflush();
}

pflush(){
   int i,j;
   for (j=0;j<256;j++){
      for (i=0;i<127;i++) {
         nodes[j][i].count[0] = 0.2;
         nodes[j][i].count[1] = 0.2;
         nodes[j][i].next[0] = &nodes[j][2*i + 1];
         nodes[j][i].next[1] = &nodes[j][2*i + 2];
      }
      for (i=127;i<255;i++) {
         nodes[j][i].count[0] = 0.2;
         nodes[j][i].count[1] = 0.2;
         nodes[j][i].next[0] = &nodes[i+1][0];
         nodes[j][i].next[1] = &nodes[i-127][0];
      }
   }
   nodesptr = nodebuf;
   preset();
}

preset(){
   p = &nodes[0][0];
}

float predict(){
   return   p->count[0] / (p->count[0] + p->count[1]);
}

pupdate(b)
   int b;
{
   float r;
   if (p->count[b] >= threshold &&
      p->next[b]->count[0]+p->next[b]->count[1]
       >= bigthresh + p->count[b]){
      new = nodesptr++;
      p->next[b]->count[0] -= new->count[0] =
         p->next[b]->count[0] * 
         (r = p->count[b]/(p->next[b]->count[1]+p->next[b]->count[0]));
      p->next[b]->count[1] -= new->count[1] =
         p->next[b]->count[1] * r;
      new->next[0] = p->next[b]->next[0];
      new->next[1] = p->next[b]->next[1];
      p->next[b] = new;
   }
   p->count[b]++;
   p = p->next[b];
   if (nodesptr > nodemaxp){
      fprintf(stderr,"flushing ...\n");
      pflush();
   }
}
