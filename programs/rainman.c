#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define EPS 1e-12
#define D 20
#define ZEROSTOP .1


int depth,nodenum,edgenum,*size,**node,*insize,**in,*outsize,**out,**inedge,**outedge;
double *w,*val,*outgrad,*bias,*bias0,*sig,*sigvel,*invec;
int vecsize,classnum,symsize,amax,pos[256];
int itercount,falsecount;
FILE *trainptr,*testptr;
fpos_t trainstart,teststart;


int opendata(char *trainfile, char *testfile)
{
int i;
char sym,sym2;
int symsize2,amax2,datasize,datasize2,classnum2;

trainptr=fopen(trainfile,"r");
if(!trainptr)
	{
	printf("trainfile not found\n");
	return 0;
	}
	
testptr=fopen(testfile,"r");
if(!testptr)
	{
	printf("testfile not found\n");
	return 0;
	}
	
fscanf(trainptr,"%d",&symsize);
fscanf(testptr,"%d",&symsize2);

if(symsize!=symsize2)
	goto err;

if(symsize==1)
	{
	fscanf(trainptr,"%d%d",&amax,&datasize);
	fscanf(testptr,"%d%d",&amax2,&datasize2);
	
	if(amax!=amax2 || datasize!=datasize2)
		goto err;
	
	vecsize=2*datasize;
	}
else
	{
	for(i=0;i<symsize;++i)
		{
		fscanf(trainptr," %c",&sym);
		fscanf(testptr," %c",&sym2);
		
		if(sym!=sym2)
			goto err;
		
		pos[sym]=i;
		}
		
	fscanf(trainptr,"%d",&datasize);
	fscanf(testptr,"%d",&datasize2);
	
	if(datasize!=datasize2)
		goto err;
	
	vecsize=symsize*datasize;
	}
	
invec=malloc(vecsize*sizeof(double));

fscanf(trainptr,"%d",&classnum);
fscanf(testptr,"%d",&classnum2);

if(classnum!=classnum2)
	goto err;

fgetpos(trainptr,&trainstart);
fgetpos(testptr,&teststart);
	
return 1;
	
err:

printf("trainfile and testfile not compatible\n");
return 0;
}


void closedata()
{
fclose(trainptr);
fclose(testptr);
}


int readdata(FILE *fp,double *vec,int *class)
{
int i,j;
char sym;
int anum;

if(symsize==1)
	for(j=0;j<vecsize/2;++j)
		{
		if(fscanf(fp," %d",&anum)==EOF)
			return 1;
			
		vec[2*j]=((double)anum)/amax;
		vec[2*j+1]=1.-vec[2*j];
		}
else
	{
	for(i=0;i<vecsize;++i)
		vec[i]=0.;

	for(j=0;j<vecsize/symsize;++j)
		{
		if(fscanf(fp," %c",&sym)==EOF)
			return 1;
			
		vec[symsize*j+pos[sym]]=1.;
		}
	}
	
fscanf(fp,"%d",class);

return 0;
}


int getnet(char *netfile)
{
int d,k,n,n1,n2,e,hidden;
FILE *fp;

fp=fopen(netfile,"r");
if(!fp)
	{
	printf("netfile not found\n");
	return 0;
	}

fscanf(fp,"%d",&hidden);
depth=hidden+1;

size=malloc((depth+1)*sizeof(int));

node=malloc((depth+1)*sizeof(int *));

nodenum=0;
for(d=0;d<=depth;++d)
	{
	fscanf(fp,"%d",&size[d]);
	
	node[d]=malloc(size[d]*sizeof(int));
	
	for(k=0;k<size[d];++k)
		node[d][k]=nodenum++;
	}

if(size[0]!=vecsize || size[depth]!=classnum)
	{
	printf("net must have %d input nodes and %d output nodes\n",vecsize,classnum);
	
	fclose(fp);
	
	return 0;
	}
	
insize=malloc(nodenum*sizeof(int));
outsize=malloc(nodenum*sizeof(int));

in=malloc(nodenum*sizeof(int*));
out=malloc(nodenum*sizeof(int*));

inedge=malloc(nodenum*sizeof(int*));
outedge=malloc(nodenum*sizeof(int*));

w=malloc(nodenum*sizeof(double));
val=malloc(nodenum*sizeof(double));
outgrad=malloc(nodenum*sizeof(double));
	
for(n=0;n<nodenum;++n)
	{
	insize[n]=0;
	outsize[n]=0;
	}

fscanf(fp,"%d",&edgenum);

bias=malloc(edgenum*sizeof(double));
bias0=malloc(edgenum*sizeof(double));
sig=malloc(edgenum*sizeof(double));
sigvel=malloc(edgenum*sizeof(double));

for(e=0;e<edgenum;++e)
	{
	fscanf(fp,"%d%d%lf",&n1,&n2,&bias[e]);
	++insize[n2];
	++outsize[n1];
	}

for(n=0;n<nodenum;++n)
	{
	in[n]=malloc(insize[n]*sizeof(int));
	out[n]=malloc(outsize[n]*sizeof(int));
	
	inedge[n]=malloc(insize[n]*sizeof(int));
	outedge[n]=malloc(outsize[n]*sizeof(int));
	}
	
fclose(fp);

fp=fopen(netfile,"r");

fscanf(fp,"%*d");

for(d=0;d<=depth;++d)
	fscanf(fp,"%*d");
	
fscanf(fp,"%*d");

for(n=0;n<nodenum;++n)
	{
	insize[n]=0;
	outsize[n]=0;
	}

for(e=0;e<edgenum;++e)
	{
	fscanf(fp,"%d%d%lf",&n1,&n2,&bias[e]);
	
	in[n2][insize[n2]]=n1;
	out[n1][outsize[n1]]=n2;
	
	inedge[n2][insize[n2]]=e;
	outedge[n1][outsize[n1]]=e;
	
	++insize[n2];
	++outsize[n1];
	}

fclose(fp);

return 1;
}


void setinput(double *input)
{
int k;

for(k=0;k<size[0];++k)
	val[node[0][k]]=input[k];
}


void setoutgrad(int class)
{
int k;

for(k=0;k<size[depth];++k)
	outgrad[node[depth][k]]=0.;
	
outgrad[node[depth][class]]=1.;
}


void calcsig(int class, int *true, int *zero, int *mult)
{
int d,k,i,outnode,innode,edge;
double s,outclass,minval;

for(d=1;d<=depth;++d)
	{
	for(k=0;k<size[d];++k)
		{
		outnode=node[d][k];
		val[outnode]=0.;
		
		for(i=0;i<insize[outnode];++i)
			{
			innode=in[outnode][i];
			edge=inedge[outnode][i];
			
			s=val[innode]-bias[edge];
			sig[edge]=s>EPS ? s : 0.;
				
			val[outnode]+=sig[edge];
			}
			
		val[outnode]*=w[outnode];
		}
	}

outclass=val[node[depth][class]];
*true=1;
*zero= outclass<EPS ? 1 : 0;
*mult=1;
	
for(k=0;k<size[depth];++k)
	{
	if(k==class)
		continue;
		
	outnode=node[depth][k];
	
	if(fabs(val[outnode]-outclass)<EPS)
		++(*mult);
	else if(val[outnode]<outclass)
		{
		*true=0;
		return;
		}
	}
		
return;
}

		
void calcoutgrad()
{
int d,k,i,innode;

for(d=depth-1;d>0;--d)
	{
	for(k=0;k<size[d];++k)
		{
		innode=node[d][k];
		
		outgrad[innode]=0.;
		
		for(i=0;i<outsize[innode];++i)
			if(sig[outedge[innode][i]]>0.)
				outgrad[innode]+=outgrad[out[innode][i]];
			
		outgrad[innode]*=w[innode];
		}
	}
}


void calcsigvel()
{
int d,k,i,j,innode,outnode,edgeout,edgein;

for(d=0;d<depth;++d)
	{
	for(k=0;k<size[d];++k)
		{
		innode=node[d][k];
	
		for(i=0;i<outsize[innode];++i)
			{
			outnode=out[innode][i];
			edgeout=outedge[innode][i];
			
			sigvel[edgeout]=0.;
			if(sig[edgeout]==0.)
				continue;
			
			if(d>0)
				{
				for(j=0;j<insize[innode];++j)
					{
					edgein=inedge[innode][j];
					if(sig[edgein]>0.)
						sigvel[edgeout]+=sigvel[edgein];
					}
					
				sigvel[edgeout]*=w[innode];
				}
			
			sigvel[edgeout]+=outgrad[outnode];
			}
		}
	}
}
			

void incbias()
{
int d,k,i,edge,outnode;
double r,rmin;

rmin=DBL_MAX;

for(edge=0;edge<edgenum;++edge)
	if(sig[edge]>0.)
		{
		r=sig[edge]/sigvel[edge];
		if(r<rmin)
			rmin=r;
		}

for(d=1;d<=depth;++d)
for(k=0;k<size[d];++k)
	{
	outnode=node[d][k];
	for(i=0;i<insize[outnode];++i)
		{
		edge=inedge[outnode][i];
		if(sig[edge]>0.)
			bias[edge]+=rmin*outgrad[outnode];
		}
	}
}


void initbias()
{
int e;

for(e=0;e<edgenum;++e)
	bias[e]=0.;
}


void setweights()
{
int d,k,n;

for(d=1;d<depth;++d)
for(k=0;k<size[d];++k)
	{
	n=node[d][k];
	w[n]=sqrt(1./(insize[n]*outsize[n]));
	}
	
for(k=0;k<size[depth];++k)
	{
	n=node[depth][k];
	w[n]=1.;
	}
}


void reducebias()
{
int d,k,i,edge,outnode,innode;

for(d=1;d<=depth;++d)
for(k=0;k<size[d];++k)
	{
	outnode=node[d][k];
	for(i=0;i<insize[outnode];++i)
		{
		edge=inedge[outnode][i];
		innode=in[outnode][i];
		
		if(bias[edge]>val[innode])
			bias[edge]=val[innode]>bias0[edge] ? val[innode] : bias0[edge];
		}
	}
}


int learn(double *input, int class)
{
int e,iter,done,true,zero,mult;

for(e=0;e<edgenum;++e)
	bias0[e]=bias[e];
	
setinput(input);
setoutgrad(class);

iter=0;

done=0;
while(!done)
	{
	calcsig(class, &true, &zero, &mult);
	
	if(zero || (true && mult==1))
		done=1;
	else
		{
		calcoutgrad();
		calcsigvel();
		incbias();
	
		++iter;
		}
	}
	
if(iter)
	reducebias();

return iter;
}


double train(int batch)
{
int t,class,iter,false;
double aveiter;

aveiter=0.;
false=0;

t=0;
while(t<batch)
	{
	if(readdata(trainptr,invec,&class))
		{
		fsetpos(trainptr,&trainstart);
		continue;
		}
	
	++t;
	
	iter=learn(invec,class);
	
	itercount+=iter;
	
	aveiter+=iter;
	false+=iter>0 ? 1 : 0;
	}
	
falsecount+=false;

return false ? aveiter/false : 0.;
}


void test(int batch, double *avetrue, double *avezero, double aveact[D])
{
int d,t,i,class,true,zero,mult,k,innode,actcount,edgecount;

*avetrue=0.;
*avezero=0.;

for(d=0;d<depth;++d)
	aveact[d]=0.;
	
fsetpos(testptr,&teststart);

t=0;
while(t<batch)
	{
	if(readdata(testptr,invec,&class))
		{
		fsetpos(testptr,&teststart);
		break;
		}
		
	++t;
	
	setinput(invec);
	
	calcsig(class, &true, &zero, &mult);
	
	*avetrue+=true ? 1./mult : 0.;
	
	if(zero && mult>1)
		*avezero+=1.;
	
	for(d=0;d<depth;++d)
		{
		actcount=0;
		edgecount=0;
		
		for(k=0;k<size[d];++k)
			{
			innode=node[d][k];
	
			if(val[innode]<EPS)
				continue;
			
			for(i=0;i<outsize[innode];++i)
				{
				actcount+=sig[outedge[innode][i]]==0. ? 1 : 0;
				++edgecount;
				}
			}
			
		aveact[d]+=((double)actcount)/edgecount;
		}
	}
	
*avetrue/=t;
*avezero/=t;

for(d=0;d<depth;++d)
	aveact[d]/=t;
}


void printnet(char *netfile)
{
FILE *fp;
int d,k,i,innode,outnode,edgeout;

fp=fopen(netfile,"w");

fprintf(fp,"%d\n",depth-1);

for(d=0;d<=depth;++d)
	fprintf(fp,"%d ",size[d]);
	
fprintf(fp,"\n%d\n\n",edgenum);

for(d=0;d<depth;++d)
for(k=0;k<size[d];++k)
	{
	innode=node[d][k];
	
	for(i=0;i<outsize[innode];++i)
		{
		outnode=out[innode][i];
		edgeout=outedge[innode][i];
			
		fprintf(fp,"%d %d %15.10e\n",innode,outnode,bias[edgeout]);
		}
	}

fclose(fp);
}


int main(int argc,char* argv[])
{
char *trainfile,*testfile,*netfile,logfile[50];
int d,t,trainstop,trainbatch,testbatch;
double avetrue,avezero,aveiter,aveact[D],avemax;
FILE *fp;

if(argc==7)
	{
	trainfile=argv[1];
	testfile=argv[2];
	netfile=argv[3];
	trainstop=atoi(argv[4]);
	trainbatch=atoi(argv[5]);
	testbatch=atoi(argv[6]);
	}
else
	{
	fprintf(stderr,"expected six arguments: trainfile, testfile, netfile, trainstop, trainbatch, testbatch\n");
	return 1;
	}
	
strcpy(logfile,netfile);
strcat(logfile,".log");

if(!opendata(trainfile,testfile))
	return 1;

if(!getnet(netfile))
	return 1;

fp=fopen(logfile,"w");

fprintf(fp,"trainfile: %s\n",trainfile);
fprintf(fp,"testfile:  %s\n\n",testfile);

fprintf(fp,"      data     false      iter aveiter   true   zero");
for(d=0;d<depth;++d)
	fprintf(fp,"    act %2d",d);
fprintf(fp,"\n");

fclose(fp);

setweights();

initbias();

itercount=0;
falsecount=0;

avemax=0.;

for(t=0;t<trainstop;t+=trainbatch)
	{
	aveiter=train(trainbatch);
	
	test(testbatch, &avetrue, &avezero, aveact);
	
	fp=fopen(logfile,"a");
	fprintf(fp,"%10.2e%10.2e%10.2e%8.2lf%7.2lf%7.2lf",(double)(t+trainbatch),(double)falsecount,(double)itercount,aveiter,100*avetrue,100*avezero);
	for(d=0;d<depth;++d)
		fprintf(fp,"%10.2e",1.-aveact[d]);
	fprintf(fp,"\n");
	fclose(fp);
	
	if(avetrue>avemax)
		{
		avemax=avetrue;
		printnet(netfile);
		}
		
	if(aveiter==0. || avezero>ZEROSTOP)
		break;
	}

fp=fopen(logfile,"a");
fprintf(fp,"\nbest:%7.2lf\ntotal false:%10d\ntotal iterations:%10d\n",100*avemax,falsecount,itercount);
fclose(fp);

closedata();

return 0;
}