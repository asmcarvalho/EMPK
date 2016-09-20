package pkcluster;

public class EMAlgorithm {
	int n;
	int N;
	int M;
	int[] cluster;
	double[] T;
	double[][] Y;
	double[][] C;
	double[][] AB;
	double[] Sigma;
	double conv;
	double[]w;
	
	public EMAlgorithm(double[] AT, double[][] AY, int M0){
		int i;
		double m;
				
		n=AT.length;
		N=AY.length;
		if(M0==0)M=N/3;
		else M=M0;
//		if(M<10)M=10;
		T=AT;
		Y=AY;
		conv=0;
		
		Sigma=new double[M];
		for(i=0;i<M;i++)Sigma[i]=1;
		
		w=new double[M];
		m=1;m=m/M;
		for(i=0;i<M;i++)w[i]=m;
		
		C=new double[M][3];
		
		//DEFINIDOS NO FIM
		cluster=new int[N];
		AB=new double[N][3];
		}
	public double integral(double[] C, double max){
		double I=0;
		double a,b1,b2;
		a=C[0];b1=C[1];b2=C[2];
		I=a*((1-Math.exp(-b1*max))/b1+(Math.exp(-b2*max)-1)/2);
		return I;
	}
	public double[][] gera(int M0){
		double[][] C=new double[M0][3]; 
		int l;
		double a,b1,b2;
		
		for(l=0;l<M0;l++){
			a=75*Math.random();
			b1=0.1+2*Math.random();
			b2=0.1+2*Math.random();
			
			C[l][0]=a;
			if(b1<b2){C[l][1]=b1;C[l][2]=b2;}
			else{C[l][1]=b2;C[l][2]=b1;}
			}
		return C;		
	}
	public static double[] gera0(){
		double[] C=new double[3]; 
		double a,b1,b2;
		
		a=75*Math.random();
		b1=0.1+2*Math.random();
		b2=0.1+2*Math.random();
		C[0]=a;
		if(b1<b2){C[1]=b1;C[2]=b2;}
		else{C[1]=b2;C[2]=b1;}
		
		return C;		
	}
	public void inicialize(Data A){
		int l,bool=0;
		double area=A.area();
		double I;
		C=gera(M);
		while(bool==0){
			bool=1;
			for(l=0;l<C.length;l++){
				I=integral(C[l],1.5*A.T[A.T.length-1]);
				if(!(area*0.2 < I && I < area*2/4)){
					bool=0;
					C[l]=gera0();
				}
			}
		}
	}
	
	public double[][] actF(){
		double[][] F=new double[M][n];
		int l,j;
		for(l=0;l<M;l++)for(j=0;j<n;j++)
			F[l][j]=C[l][0]*(Math.exp(-C[l][1]*T[j])-Math.exp(-C[l][2]*T[j]));
	return F;
	}
	public double Pil(double[][]F,int i, int l){
		double r,s,y,x,p;
		int j;
		
		r=Math.PI*Sigma[l]*2;
		r=Math.log(r);
		r=-r*0.5*n;
		s=0;
		for(j=0;j<n;j++){
			y=Y[i][j]-F[l][j];
			y=y*y;
			s=s+y;}
		x=s/(Sigma[l]*2);
		p=r-x;
		
		return p;
	}
	public double[][] actP(double[][] F){
		//CALCULADO LOG(P)
		
		double[][] P=new double[N][M];
		int l,i;
		double p;
		
		for(i=0;i<N;i++)
			for(l=0;l<M;l++){
				p=Pil(F,i,l);
				P[i][l]=p;
			}
		return P;
	}
	public double[] actW(double[][]P){
		//CALCULADO LOG(W)
		double[]W=new double[N];
		int i,l;
		double x,y;
		for(i=0;i<N;i++){
			y=0;
			for(l=0;l<M;l++){
				x=Math.log(w[l])+P[i][l];
				//PARA NÃO DAR 0
				x=x+500;
				x=Math.exp(x);
				y=y+x;
			}
			W[i]=Math.log(y);

		}
		return W;
	}
	public double[][] actX(double[][]P,double[]W){
		//CALCULADO LOG(X)
		//GUARDADO X
		double[][]X=new double[N][M];
		int i,l;
		double ps=Math.pow(10, -10);
		double x;
		double[] Xt=new double[M];
		for(l=0;l<M;l++)Xt[l]=0;
		
		for(i=0;i<N;i++)for(l=0;l<M;l++){
				//500 para evitar os NaN em EXP
				x=Math.log(w[l])+P[i][l]+500;
				x=Math.exp(x);
				//ps para evitar uma coluna inteira como 0
				x=x+ps;
				x=x/(Math.exp(W[i])+M*ps);
				
				X[i][l]=x;
				Xt[l]=Xt[l]+x;}
		
		//ACTUALIZA PESOS w
		w=new double[M];

		for(l=0;l<M;l++){
			w[l]=Xt[l];
			w[l]=w[l]/N;}
		
		return X;
	}
	public void actSigma(double[][]X,double[][]F){
		double x,y;
		int i,j,l;
		for(l=0;l<M;l++){
			x=0;
			y=0;
			for(i=0;i<N;i++){
				y=y+n*X[i][l];
				for(j=0;j<n;j++){
					x=x+X[i][l]*(Y[i][j]-F[l][j])*(Y[i][j]-F[l][j]);
				}
			}
			
		Sigma[l]=x/y;
		}
	}
	//if(cov==1){FIM}else{continua}
	public double zeroa(int l,double[][]X,double[][]F){
		double a,x,y,z;
		int i,j;
		x=0;
		y=0;
		for(j=0;j<n;j++){
			z=Math.exp(-C[l][1]*T[j])-Math.exp(-C[l][2]*T[j]);
			for(i=0;i<N;i++){
				x=x+X[i][l]*z*Y[i][j];
				y=y+X[i][l]*z*z;
				}}
		a=x/y;
		return a;
	}
	public void actC(double[][]X,double[][]F){
		double a,b1,b2;
		double x1;
		double[] bn=new double[2];
		int l,bool;
		bool=1;
		for(l=0;l<M;l++){
			a=zeroa(l,X,F);C[l][0]=a;
			//STOP B1
			bn=actB(C[l][1],C[l][2],l,X);
			b1=bn[0];
			b2=bn[1];
			
			x1=b1-C[l][1];
			x1=x1*x1;
			if(x1>0.000001*1)bool=0;
			//STOP B2
			x1=b2-C[l][2];
			x1=x1*x1;
			if(x1>0.000001*1)bool=0;
			//UPDATE
			C[l][1]=b1;
			C[l][2]=b2;
			}
		if(bool==1)conv=1;
	}
	public int maximo(int i,double[][]X){
		int id=-1;
		double m=0;
		int l;
		for(l=0;l<M;l++)if(X[i][l]>m){m=X[i][l];id=l;}
		
/*ASMC		if(id==-1){
			System.out.print("ID MAXIMO: "+M);
			for(l=0;l<M;l++)System.out.print(X[i][l]+"; ");
			System.out.println();
			for(l=0;l<M;l++)System.out.println(C[l][0]+"; "+C[l][1]+"; "+C[l][2]);
		}
*/		
		return id;
	}
	public int[] actualize(double[][]X){
		int i,l;
		int[] cmax=new int[M];
		for(l=0;l<M;l++)cmax[l]=0;
		for(i=0;i<N;i++){
			l=maximo(i,X);
			if(l!=-1){
				AB[i][0]=C[l][0];
				AB[i][1]=C[l][1];
				AB[i][2]=C[l][2];
				cluster[i]=l;}
			else{ cluster[i]=0;}
			cmax[l]++;
						
		}
		return cmax;
	}
	
	public int colapsa(double D, double omega, int[] cmax){
		int bool=0;
		double x1,x2,x,d;
		int j,l,l1,m;
		d=D;
		m=0;
		//COLAPSO 1
		for(l=0;l<M;l++)for(l1=l+1;l1<M;l1++)if(w[l]*w[l1]!=0){
			x=0;
			for(j=0;j<n;j++){
				x1=C[l][0]*(Math.exp(-C[l][1]*T[j])-Math.exp(-C[l][2]*T[j]));
				x2=C[l1][0]*(Math.exp(-C[l1][1]*T[j])-Math.exp(-C[l1][2]*T[j]));
				x=x+(x1-x2)*(x1-x2);
			}
			if(x<d){
				C[l][0]=(C[l][0]+C[l1][0])/2;
				C[l][1]=(C[l][1]+C[l1][1])/2;
				C[l][2]=(C[l][2]+C[l1][2])/2;
				w[l]=w[l]+w[l1];
				Sigma[l]=Sigma[l]+Sigma[l1];
				w[l1]=0;
				cmax[l]=cmax[l]+cmax[l1];
				m++;
			}
		}
		
		//COLAPSO 2
		for(l=0;l<M;l++)if(w[l]!=0&&w[l]<omega
				//asmc&&cmax[l]==0
				){
			w[l]=0;
			m++;
		}
		
		//ACTUALIZA OS PARAMETORS APOS COLAPSAR, SE COLAPSOU
		if(m>0){
//			System.out.println("COLAPSA");
			conv=0;
			bool=1;
			double[] w1=new double[M-m];
			double[][] C1=new double[M-m][3];
			double[] s=new double[M-m];
			j=0;
			for(l=0;l<M;l++)if(w[l]!=0){
				w1[j]=w[l];
				C1[j][0]=C[l][0];
				C1[j][1]=C[l][1];
				C1[j][2]=C[l][2];
				s[j]=Sigma[l];
				j++;
			}
			C=C1;
			M=M-m;
			w=w1;
			Sigma=s;
		}
		//SE NAO ALTERA, bool=0 => PARA O ALGORITMO
		return bool;
	}
	public double likel(double[][] P0, double[][] X0){
		double res=0;
		double[][] F;
		double[][] P;
		double[] W;
		double[][] X;
		int i,l;
		
		F=actF();
		P=actP(F);
		W=actW(P);
		X=actX(P,W);
		
		for(i=0;i<N;i++)for(l=0;l<M;l++)
			res=res+X[i][l]*(Math.log(w[l])+P[i][l]);
				
		return res;
	}
	
	public double[] Newton(double b1, double b2, int l, double[][] X){
		double a=C[l][0];
		double h1=0,h2=0,dh1=0,dh2=0,ddh1=0,ddh2=0;
		double[] res= new double[6];
		
		int i,j;
		double e1,e2,aux,bux,k0,k1,k2;
		for(j=0;j<n;j++){
			e1=Math.exp(-b1*T[j]);
			e2=Math.exp(-b2*T[j]);
			bux=a*(e1-e2);
			for(i=0;i<N;i++){
				k0=X[i][l]*T[j];
				k1=k0*T[j];
				k2=k1*T[j];
				aux=Y[i][j]-bux;
				h1=h1+aux*k0*e1;
				h2=h2+aux*k0*e2;
				dh1=dh1+k1*e1*(a*e1-aux);
				dh2=dh2+k1*e2*(-a*e2-aux);
				ddh1=ddh1+k2*e1*(3*a*e1-aux);
				ddh2=ddh2+k2*e2*(-3*a*e2-aux);
			}
		}
		
		h1=-h1*a/Sigma[l];
		dh1=-dh1*a/Sigma[l];
		h2=h2*a/Sigma[l];
		dh2=dh2*a/Sigma[l];
		ddh1=ddh1*a/Sigma[l];
		ddh2=-ddh2*a/Sigma[l];
		
		res[0]=h1;res[1]=dh1;res[2]=ddh1;
		res[3]=h2;res[4]=dh2;res[5]=ddh2;
		return res;
	}
	public double[] Itera(double b1, double b2, int l, double[][] X, int bo1, int bo2, int max,double par, double tol1, double tol2){
		double[] res=new double[4];
		
		double[] bn= new double[2];
		double[] b= new double[2];
		double[] L;
		int bool1,bool2,k;
		double x1,x2,h1,h2;
		k=0;
		bool1=bo1;
		bool2=bo1;
		bn[0]=b1;
		bn[1]=b2;
		
		while(k<max&&(bool1==1||bool2==1)){
			b[0]=bn[0];b[1]=bn[1];
			L=Newton(b[0],b[1],l,X);
			h1=L[0];h2=L[3];
			if(bool1==1&&L[1]!=0)bn[0]=b[0]-h1/L[1]*par;
			if(bool2==1&&L[4]!=0)bn[1]=b[1]-h2/L[4]*par;
			x1=b[0]-bn[0];
			x2=b[1]-bn[1];
			if(x1<0)x1=-x1;
			if(x2<0)x2=-x2;
			if(h1<0)h1=-h1;
			if(h2<0)h2=-h2;
			
			if(x1<tol1||h1<tol2)bool1=0;
			if(x2<tol1||h2<tol2)bool2=0;
			
			k++;
		}
		
		res[0]=bn[0];res[1]=bn[1];res[2]=bool1;res[3]=bool2;
		
		return res;
	}
	public double[] Knee(double b1, double b2, int l, double[][] X, int bo1, int bo2, int max,double par, double tol1, double tol2){
		double[] res=new double[4];
		
		double[] bn= new double[2];
		double[] b= new double[2];
		double[] L;
		int bool1,bool2,k;
		double x1,x2,h1,h2,d1,d2;
		k=0;
		bool1=bo1;
		bool2=bo1;
		bn[0]=b1;
		bn[1]=b2;
		
		while(k<max&&(bool1==1||bool2==1)){
			b[0]=bn[0];b[1]=bn[1];
			L=Newton(b[0],b[1],l,X);
			h1=L[0];h2=L[3];
			d1=1+h1*L[2]+L[1]*L[1];
			d2=1+h2*L[5]+L[4]*L[4];
			if(bool1==1&&d1!=0)bn[0]=b[0]-(b[0]+h1*L[1])/d1*par;
			if(bool2==1&&d2!=0)bn[1]=b[1]-(b[1]+h2*L[4])/d2*par;
			x1=b[0]-bn[0];
			x2=b[1]-bn[1];
			if(x1<0)x1=-x1;
			if(x2<0)x2=-x2;
			if(h1<0)h1=-h1;
			if(h2<0)h2=-h2;
			
			if(x1<tol1||h1<tol2)bool1=0;
			if(x2<tol1||h2<tol2)bool2=0;
			
			k++;
		}
		
		res[0]=bn[0];res[1]=bn[1];res[2]=bool1;res[3]=bool2;
		
		return res;
	}
	public double[] actB(double b1, double b2, int l, double[][] X){
		double[] bn= new double[2];
		double[] it;
		double[] L;
		int bool1,bool2;
		double x1;
		
		it=Itera(b1,b2,l,X,1,1,100000,0.2,0.0000001,0.0000000001);
		bn[0]=it[0];bn[1]=it[1];
		bool1=(int)it[2];
		bool2=(int)it[3];
				
		if(bool1==1||bool2==1){
			L=Newton(bn[0],bn[1],l,X);
			if(L[0]<0)L[0]=-L[0];
			if(L[3]<0)L[3]=-L[3];
			if(L[0]>0.01||L[3]>0.01||L[1]>0||L[4]>0){
				it=Itera(bn[0],bn[1],l,X,bool1,bool2,10000,0.5,0.00001,0.00000001);
				bn[0]=it[0];bn[1]=it[1];
				bool1=(int)it[2];
				bool2=(int)it[3];}
			if(bool1==1||bool2==1){
				L=Newton(bn[0],bn[1],l,X);
				if(L[0]<0)L[0]=-L[0];
				if(L[3]<0)L[3]=-L[3];
//ASMC				if(L[0]>0.01||L[3]>0.01||L[1]>0||L[4]>0){
//ASMC				System.out.print("#### AUMENTAR K #### "+bool1+";"+bool2+";");
//ASMC				System.out.println(bn[0]+";"+bn[1]+";"+L[0]+";"+L[3]+";"+L[1]+";"+L[4]);}
				}
		}
		
		if(!(0<bn[0]&&bn[0]<5))bool1=2;
		if(!(0<bn[1]&&bn[1]<5))bool2=2;
		if(bn[0]>bn[1]){x1=bn[0];bn[0]=bn[1];bn[1]=x1;}
				
		if(bool1==2||bool2==2){
			if(bool1==2&&bool2==2)
				it=Knee(b1,b2,l,X,1,1,10000,0.5,0.000001,0.00000001);
			if(bool1!=2&&bool2==2)
				it=Knee(bn[0],b2,l,X,0,1,10000,0.5,0.000001,0.00000001);
			if(bool1==2&&bool2!=2)
				it=Knee(b1,bn[1],l,X,1,0,10000,0.5,0.000001,0.00000001);
			bn[0]=it[0];bn[1]=it[1];
			
			if(!(0<bn[0]&&bn[0]<5))bn[0]=b1;
			if(!(0<bn[1]&&bn[1]<5))bn[1]=b2;
			
			
//			bn=Cotovelo(b1,b2,l,X,bool1,bool2,bn);
			
			L=Newton(bn[0],bn[1],l,X);
//ASMC			System.out.print("#### COTOVELO #### ");
//ASMC			System.out.println(bn[0]+";"+bn[1]+";"+L[0]+";"+L[3]+";"+L[1]+";"+L[4]);
//ASMC			int j;
//ASMC			for(j=0;j<C.length;j++)	System.out.print("{"+C[j][0]+","+C[j][1]+","+C[j][2]+"},");
//ASMC			System.out.println();
		}
			
		if(bn[0]>bn[1]){x1=bn[0];bn[0]=bn[1];bn[1]=x1;}
		return bn;
	}
		
	public double[] Cotovelo(double b1, double b2, int l, double[][] X, int bool1, int bool2,double[]BN){
		double[] bn= new double[2];
		double[] b= new double[2];
		int k;
		double B=-0.3;
		double x1,x2;
		double[] aux=new double[4];
		double[] dh=new double[2];
		if(bool1==2)bool1=1;
		if(bool2==2)bool2=1;
		bn[0]=BN[0];
		bn[1]=BN[1];
		
		while(B>-5&&(bool1==1||bool2==1)){
			k=0;
			if(bool1==1)bn[0]=b1;
			if(bool2==1)bn[1]=b2;
			while(k<10000&&(bool1==1||bool2==1)){
				b[0]=bn[0];b[1]=bn[1];
				aux=NewtonC(b[0],b[1],l,X,bool1,bool2);
				
				bn[0]=aux[0];bn[1]=aux[1];
				dh[0]=aux[2];dh[1]=aux[3];
				
				x1=b[0]-bn[0];
				x2=b[1]-bn[1];
				if(x1<0)x1=-x1;
				if(x2<0)x2=-x2;
				if(x1<0.000001*0.0001||dh[0]*bool1>B)bool1=0;
				if(x2<0.000001*0.0001||dh[1]*bool2>B)bool2=0;
				
				if(!(0<bn[0]&&bn[0]<bn[1]))bool1=2;
				if(!(bn[0]<bn[1]&&bn[1]<5))bool1=2;
				
				k++;
			}
			if(bool1==2||bool2==2)B=B-0.2;
			if(bool1==2)bool1=1;
			if(bool2==2)bool2=1;
		}
		
		if(B<=-5){
			if(0>bn[0])bn[0]=b1;
			if(5<bn[1])bn[1]=b2;
			if(bn[0]>bn[1]){x1=bn[0];bn[0]=bn[1];bn[1]=x1;}
//			if(bool1==1)bn[0]=b1;
//			if(bool2==1)bn[1]=b2;
			}
		
//		System.out.println("#### COTOVELO ###");
//		if(!(0<bn[0]&&bn[0]<bn[1]&&bn[1]<5))
//			System.out.println("ERRO: "+ bn[0]+"; "+bn[1]+"; "+l+"; "+C[l][0]+"; "+B);
		
		return bn;
	}
		
	public double[] NewtonC(double b1, double b2, int l, double[][] X, int bool1, int bool2){
		double a=C[l][0];
		double h1=0,h2=0,dh1=0,dh2=0;
		double[] res= new double[4];
		res[0]=b1; res[1]=b2;
		int i,j;
		double e1,e2,aux,bux,k0,k1;
		for(j=0;j<n;j++){
			e1=Math.exp(-b1*T[j]);
			e2=Math.exp(-b2*T[j]);
			bux=a*(e1-e2);
			for(i=0;i<N;i++){
				k0=X[i][l]*T[j];
				k1=k0*T[j];
				aux=Y[i][j]-bux;
				h1=h1+aux*k0*e1;
				h2=h2+aux*k0*e2;
				dh1=dh1+k1*e1*(a*e1-aux);
				dh2=dh2+k1*e2*(-a*e2-aux);
			}
		}
		
		if(bool1==1)res[0]=res[0]-h1/dh1;
		if(bool2==1)res[1]=res[1]-h2/dh2;
		res[2]=dh1;res[3]=dh2;
		
		return res;
	}
	
	public Output runEM(){
		Output out;
		
		double D1,omega;
		//POR ANALISAR ESTES VALORES. FUTURO VEM DE INPUT
		//ASMC omega=0.05/2;
		omega=0.05/2;
		D1=9;
		
		
		int j=0,i;
		double[][] F;
		double[][] P;
		double[] W;
		double[][] X;
		int[] cmax;
		int bool=1;
		double Q=-150;
		//DEPOIS REMOVER OS PRINTS INTERNOS E O i E j DOS whiles		
		i=0;
		while(bool==1&&i</*ASMC 10*/50){
			F=actF();
			P=actP(F);
			W=actW(P);
			X=actX(P,W);
			actSigma(X,F);
			
			j=0;
			while(conv==0&&j<500){
//				System.out.println("Passo "+j);
				j++;
				actC(X,F);
				F=actF();
				P=actP(F);
				W=actW(P);
				X=actX(P,W);
				actSigma(X,F);
			}
			
			
//			System.out.println("SAI "+j);
			cmax=actualize(X);
			bool=colapsa(D1,omega,cmax);
			i++;
			if(bool==0)Q=likel(P,X);
		}
		
		
		out=new Output(C,Q,j,Sigma,w,M,cluster);
		
		
//		System.out.println("OUTPUT");
//		System.out.println();
//		
//		System.out.println("Sigma: ");
//		for(l=0;l<M;l++)System.out.print(Sigma[l]+"; ");
//		System.out.println();
//		
//		System.out.println("Pesos: ");
//		for(l=0;l<M;l++)System.out.print(w[l]+"; ");
//		System.out.println();
//		
//		System.out.println("Number of clusters: "+M);
//		
//		System.out.println("a; ke; ka");
//		for(l=0;l<M;l++)System.out.println(C[l][0]+"; "+C[l][1]+"; "+C[l][2]);
		
//		System.out.println("Dados dos pacientes");
//		System.out.println("cl; a; ke; ka");
//		for(i=0;i<N;i++)System.out.println(cluster[i]+"; "+100*AB[i][0]+"; "+AB[i][1]+"; "+AB[i][2]);
		
//		System.out.println("Q-value: "+likel());
		
		return out;
	}

	
}

