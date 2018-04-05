#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(int argc,char* argv[])
{
FILE *fp;
char *netfile;
int size[21],insize,outsize,hiddenlayers,growth,outdeg,d,depth,edgenum;
int *notused,instart,outstart,i,j,count,in;

srand(time(0));

if(argc==5 || argc==6)
	{
	insize=atoi(argv[1]);
	outsize=atoi(argv[2]);
	hiddenlayers=atoi(argv[3]);
	growth=atoi(argv[4]);
	}
else
	{
	fprintf(stderr,"expected four or five arguments: insize outsize hiddenlayers growth [netfile]\n");
	return 1;
	}
	
outdeg=2*growth;
	
if(hiddenlayers>20)
	{
	printf("hiddenlayers may not exceed 20\n");
	return 0;
	}
	
if(argc==6)
	netfile=argv[5];

depth=hiddenlayers+1;

size[0]=insize;
size[depth]=outsize;
edgenum=0;

for(d=1;d<depth;++d)
	{
	size[d]=(outdeg/2)*size[d-1];
	edgenum+=2*size[d];
	}
	
edgenum+=outsize*size[depth-1];

if(argc==5)
	{
	printf("%d edges\n",edgenum);
	return 0;
	}

notused=malloc(size[depth-1]*sizeof(int));

fp=fopen(netfile,"w");

fprintf(fp,"%d\n",hiddenlayers);

for(d=0;d<=depth;++d)
	fprintf(fp,"%d ",size[d]);
	
fprintf(fp,"\n%d\n\n",edgenum);

instart=0;
outstart=size[0];

count=0;

for(d=1;d<depth;++d)
	{
	for(in=0;in<2;++in)
	for(j=0;j<size[d];++j)
		{
		if(count==0)
			{
			for(i=0;i<size[d-1];++i)
				notused[i]=i;
		
			count=size[d-1];
			}
		
		i=rand()%count;
	
		fprintf(fp,"%d %d 0.\n",instart+notused[i],outstart+j);
		
		count--;
		notused[i]=notused[count];
		}
		
	instart+=size[d-1];
	outstart+=size[d];
	}
	
for(j=0;j<size[depth];++j)
for(i=0;i<size[depth-1];++i)
	fprintf(fp,"%d %d 0.\n",instart+i,outstart+j);
	
fclose(fp);

return 0;
}
