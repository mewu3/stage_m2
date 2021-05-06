#include "select.h"
#include "visible.h"

static void localverbosealphabet(Alphabet *alpha)
{
  Uchar c;
  Uint i;

  printf("alphabet of size %lu: ",(Showuint) alpha->mapsize);
  for(i=0; i< alpha->mapsize; i++)
  {
    c = alpha->characters[i];
    if(INVISIBLE(c))
    {
      printf("\\%lu",(Showuint) c);
    } else
    {
      printf("%c",c);
    }
  }
  printf("\n");
}

static void localshowmultiseq(Multiseq *multiseq)
{
  printf("numofsequences=%lu\n",(Showuint) multiseq->numofsequences);
  printf("totallength=%lu\n",(Showuint) multiseq->totallength);
  printf("totalnumoffiles=%lu\n",(Showuint) multiseq->totalnumoffiles);
}

static void localshowmatch(StoreMatch *storematch)
{
  printf("idnumber=%lu\n",(Showuint) storematch->idnumber);
  printf("distance=%ld\n",(Showsint) storematch->Storedistance);
  printf("pos1=%lu\n",(Showuint) storematch->Storeposition1);
  printf("len1=%lu\n",(Showuint) storematch->Storelength1);
  printf("pos2=%lu\n",(Showuint) storematch->Storeposition2);
  printf("len2=%lu\n",(Showuint) storematch->Storelength2);
  printf("seqnum1=%lu\n",(Showuint) storematch->Storeseqnum1);
  printf("relpos1=%lu\n",(Showuint) storematch->Storerelpos1);
  printf("seqnum2=%lu\n",(Showuint) storematch->Storeseqnum2);
  printf("relpos2=%lu\n",(Showuint) storematch->Storerelpos2);
  printf("evalue=%2e\n",storematch->StoreEvalue);
}

Sint selectmatch(Alphabet *alpha,
                 Multiseq *virtualmultiseq,
                 Multiseq *querymultiseq,
                 StoreMatch *storematch)
{
  printf("selectnone\n");
  localverbosealphabet(alpha);
  localshowmultiseq(virtualmultiseq);
  localshowmatch(storematch);
  return 0;
}
