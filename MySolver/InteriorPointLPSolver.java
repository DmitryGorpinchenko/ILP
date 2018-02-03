import java.io.*;

/*********************************************************************************
 *  class for solving a Linear Programs by a central path interior point method
 *  
 *  compilation: javac InteriorPointLPSolver.java  
 *         
 *  usage:       java InteriorPointLPSolver filename -t trace_flag       
 *  
 *  Author:      Gorpinchenko Dmitry, KINR NASU, 2014, Copyright
 *
 *  assumes that LP has a standard form: 
 *
 *  maximize    c*x 
 *  subject to  Ax <= b
 *  		    x >= 0
 *
 *  LP provided with a file 'filename' written in the dictionary format 
 *
 *  m n 
 *  B1 B2 ..... Bm  [indices for slack variables] 
 *  N1 N2 ..... Nn  [indices for decision variables] 
 *  b1 b2 ..... bm  [b vector components] 
 *  -a11 .... -a1n  [first row of -A matrix] 
 *       .... 
 *  -am1 .... -amn  [mth row of -A matrix] 
 *  z0 c1 ..... cn  [objective coefficients (z0 = 0.0 in the original dictionary)] 
 *
 *  Note: in brackets some comments are given, there is no such things in the file
 *********************************************************************************/
 
class InteriorPointLPSolver {
	
	public static final int    MAXITS = 200;
	public static final double EPS    = 1e-6;
	public static final double DELTA  = 0.02;
	public static final double SIGMA  = 0.9;
	public static final double BIG    = Double.POSITIVE_INFINITY;
	public static final double BOUND  = 1e9;
	
	private double[][] A;
	private double[] b;
	private double[] c;
	private double[] x;
	private double primal_obj;
	private int m;
	private int n;	
	private int iter;
	private boolean is_optimal, is_infeasible, is_unbounded; 
	
	public InteriorPointLPSolver(double[][] A, double[] b, double[] c) {
		this.A = A;
		this.b = b;
		this.c = c;
		n = c.length;
		m = b.length;
		x = new double[n];
	}
	
	public void solve(boolean trace_flag) {
		if(trace_flag) {
			StdOut.println("\niter |  primal obj |    dual obj |      normrp |      normrd |     normgap |\n");
		}
		double[] xs  = new double[m];
		double[] y   = new double[m];
		double[] ys  = new double[n];
		double[] dx  = new double[n];
		double[] dxs = new double[m];
		double[] dy  = new double[m];
		double[] dys = new double[n];
		double[] rp  = new double[m];
		double[] rd  = new double[n];
		double[] rhs = new double[m];
		double[] ATy = new double[n];
		double[] Ax  = new double[m]; 
		double[] dn  = new double[n];
		double[] dm  = new double[m];
		double[] tempn = new double[n];
		double[] tempm = new double[m];
		double[][] tempmm = new double[m][m];
		double[][] AT = transpose(A);
		double rpfact = 1 + Math.sqrt(dotprod(b, b));
		double rdfact = 1 + Math.sqrt(dotprod(c, c));
		set_initial_point(x, xs, y, ys);
		double normrp_old = BIG;
		double normrd_old = BIG;
		double dual_obj, normrp, normrd, normgap, gap, mu, alpha_p, alpha_d;
		for(iter = 0; iter < MAXITS; iter++) {
			matvecprod(Ax, A, x);
			for(int i = 0; i < m; i++) {
				rp[i] = Ax[i] + xs[i] - b[i];
			}
			normrp = Math.sqrt(dotprod(rp, rp))/rpfact;
			matvecprod(ATy, AT, y);
			for(int i = 0; i < n; i++) {
				rd[i] = ATy[i] - ys[i] - c[i];
			}
			normrd = Math.sqrt(dotprod(rd, rd))/rdfact;
			gap = dotprod(x, ys) + dotprod(xs, y);
			mu = DELTA*gap/(n+m);
			primal_obj = dotprod(c, x);
			dual_obj = dotprod(b, y);
			normgap = gap/(1+Math.abs(primal_obj));
			if(trace_flag) {
				print_trace_info(dual_obj, normrp, normrd, normgap);
			}
			if(normrp < EPS && normrd < EPS && normgap < EPS) {
				is_optimal = true;
				return;
			}	
			if(/* normrp > 1000*normrp_old && normrp > EPS */ inf_norm(y) > BOUND) {
				is_infeasible = true;
				return;
			}
			if(/* normrd > 1000*normrd_old && normrd > EPS */ inf_norm(x) > BOUND) {
				is_unbounded = true;
				return;
			}
			for(int i = 0; i < n; i++) {
				dn[i] = x[i]/ys[i];
			}
			for(int i = 0; i < m; i++) {
				dm[i] = xs[i]/y[i];
			}
			for(int i = 0; i < n; i++) {
				tempn[i] = x[i] - mu/ys[i] + dn[i]*rd[i];
			}
			matvecprod(tempm, A, tempn);
			for(int i = 0; i < m; i++) {
				rhs[i] = rp[i] + mu/y[i] - xs[i] - tempm[i];
			}
			for(int i = 0; i < m; i++) {
				for(int j = 0; j < m; j++) {
					tempmm[i][j] = 0;
					for(int k = 0; k < n; k++) {
						tempmm[i][j] += A[i][k]*A[j][k]*dn[k];
					}
				}
			}
			for(int i = 0; i < m; i++) {
				tempmm[i][i] += dm[i];
			} 
			LDL(dy, tempmm, rhs);
			matvecprod(tempn, AT, dy);
			for(int i = 0; i < n; i++) {
				dys[i] = tempn[i] + rd[i];
			}
			for(int i = 0; i < n; i++) {
				dx[i] = -dn[i]*dys[i] + mu/ys[i] - x[i];
			}
			for(int i = 0; i < m; i++) {
				dxs[i] = -dm[i]*dy[i] + mu/y[i] - xs[i];
			} 
			alpha_p = 1;
			for(int i = 0; i < n; i++) {
				if(x[i] + alpha_p*dx[i] < 0) {
					alpha_p = -x[i]/dx[i];
				}
			}	
			for(int i = 0; i < m; i++) {
				if(xs[i] + alpha_p*dxs[i] < 0) {
					alpha_p = -xs[i]/dxs[i];
				}
			} 
			alpha_d = 1;
			for(int i = 0; i < m; i++) {
				if(y[i] + alpha_d*dy[i] < 0) {
					alpha_d = -y[i]/dy[i];
				}
			}	
			for(int i = 0; i < n; i++) {
				if(ys[i] + alpha_d*dys[i] < 0) {
					alpha_d = -ys[i]/dys[i];
				}
			}	 
			alpha_p *= SIGMA;
			alpha_d *= SIGMA;
			for(int i = 0; i < n; i++) {
				x[i]  += alpha_p*dx[i];
				ys[i] += alpha_d*dys[i];
			}
			for(int i = 0; i < m; i++) {
				xs[i] += alpha_p*dxs[i];
				y[i]  += alpha_d*dy[i];
			}
			normrp_old = normrp;
			normrd_old = normrd;
		}
	}
	
	public void set_initial_point(double[] x, double[] xs, double[] y, double[] ys) {
		for(int i = 0; i < n; i++) {
			x[i]  = 1000.0;
			ys[i] = 1000.0;
		}
		for(int i = 0; i < m; i++) {
			xs[i] = 1000.0;
			y[i]  = 1000.0;
		}
	}
	
	public void print_trace_info(double dual_obj, double normrp, double normrd, double normgap) {
		StdOut.println(String.format(java.util.Locale.UK, "%4d | %11.4e | %11.4e | %11.4e | %11.4e | %11.4e |", iter, primal_obj, dual_obj, normrp, normrd, normgap));
	}
	
	public void print_results() {
		StdOut.println("\n******* Results *******\n");
		if(is_optimal) {
			StdOut.println("Optimal solution found after " + iter + " iteration of a central path algorithm\n");
			for(int i = 0; i < n; i++) {
				if(x[i] > EPS) {
					StdOut.println(String.format(java.util.Locale.UK, " x%-3d = %11.4e", i+1, x[i]));
				}
			}
			StdOut.println(String.format(java.util.Locale.UK, "\nOptimal objective value: %.4e\n", primal_obj));
		} else if(is_infeasible) {
			StdOut.println("LP is INFEASIBLE!!!\n");
		} else if(is_unbounded) {
			StdOut.println("LP is UNBOUNDED!!!\n");
		} else {
			StdOut.println("Number of iterations exceeded the limit of " + MAXITS + "\n");
		}
	}
	
	public static double min(double a, double b) {
		return a < b ? a : b;
	}
	
	public static void LDL(double[] x, double[][] A, double[] b) {
		int n = A.length;
		double sum = 0;
		//factorization
		double[][] L = new double[n][n];
		double[] D = new double[n];
		for(int i = 0; i < n; i++) {
			L[i][i] = 1;
		}
		for(int j = 0; j < n; j++) {
			sum = 0;
			for(int k = 0; k < j; k++) {
				sum += L[j][k]*L[j][k]*D[k];
			}
			D[j] = A[j][j] - sum;
			/* if(D[j] < 1e-40) {
				D[j] = 1e128;
			} */
			for(int i = j+1; i < n; i++) {
				sum = 0;
				for(int k = 0; k < j; k++) {
					sum += L[i][k]*L[j][k]*D[k];
				}
				L[i][j] = (A[i][j] - sum)/D[j];
			}
		}
		double[] y = new double[n];
		for(int i = 0; i < n; i++) {
			sum = 0;
			for(int k = 0; k < i; k++) {
				sum += L[i][k]*y[k];
			}
			y[i] = (b[i] - sum)/L[i][i];
		}
		//backsubstitution
		for(int i = n-1; i >= 0; i--) {
			sum = 0;
			for(int k = i+1; k < n; k++) {
				sum += D[i]*L[k][i]*x[k];
			}
			x[i] = (y[i] - sum)/(L[i][i]*D[i]);
		}
	}
	
	public static double inf_norm(double[] x) {
		double norm = -BIG;
		for(int i = 0; i < x.length; i++) {
			if(norm < Math.abs(x[i])) {
				norm = Math.abs(x[i]);
			}
		}
		return norm;
	}	
	
	public static double dotprod(double[] a, double[] b) {
		double dot = 0;
		for(int i = 0; i < a.length; i++) {
			dot += a[i]*b[i];
		}
		return dot;
	}
	
	public static void matvecprod(double[] b, double[][] A, double[] x) {
		int n = A.length, m = x.length;
		for(int i = 0; i < n; i++) {
			b[i] = dotprod(A[i], x);
		}
 	}
	
	public static double[][] transpose(double[][] A) {
		int m = A.length, n = A[0].length;
		double[][] AT = new double[n][m];
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				AT[i][j] = A[j][i];
			}
		}
		return AT;
	}
	
	public static void main(String[] args) throws IOException {
		Stopwatch sw = new Stopwatch();
		BufferedReader input = new BufferedReader(new FileReader(args[0]));
		String[] tokens = input.readLine().split("\\s+");
		int m = Integer.parseInt(tokens[0]), n = Integer.parseInt(tokens[1]);
		input.readLine();input.readLine();
		double[][] A = new double[m][n]; 
		double[] b = new double[m];
		double[] c = new double[n];		
		tokens = input.readLine().split("\\s+");
		for(int i = 0; i < m; i++) {
			b[i] = Double.parseDouble(tokens[i]);
		}
		for(int i = 0; i < m; i++) {
			tokens = input.readLine().split("\\s+");
			for(int j = 0; j < n; j++) {
				A[i][j] = -Double.parseDouble(tokens[j]);
			}
		}
		tokens = input.readLine().split("\\s+");
		for(int i = 0; i < n; i++) {
			c[i] = Double.parseDouble(tokens[i+1]);
		}
		input.close();
		InteriorPointLPSolver solver = new InteriorPointLPSolver(A, b, c);
		solver.solve(Boolean.parseBoolean(args[2]));
		solver.print_results();
		StdOut.println("\nTiming results: " + sw.elapsedTime()); 
	}
}	