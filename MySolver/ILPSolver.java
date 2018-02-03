import java.io.*;

class ILPSolver {
	
	private Simplex simplex;
	private int iter;
	private static double tol = Dictionary.tol;
	
	public ILPSolver(Dictionary dict, boolean trace_flag) {
		//solving initial LP relaxation
		try {
			if(trace_flag) {
				StdOut.println("\nSolving original LP relaxation ...\n");
			}
			simplex = new Simplex(dict, false);
		} catch (Exception e) {
			String message = e.getMessage();
			if(message.equals("Linear program is INFEASIBLE!")) {
				throw new ArithmeticException("ILP is INFEASIBLE!");
			} else if(message.equals("Linear program is UNBOUNDED!")) {
				throw new ArithmeticException("ILP is UNBOUNDED!");
			}
		}
		//solving futher by adding the Gomory cuts
		if(trace_flag) {
			StdOut.println("\nSolving ILP by using Gomory cuts ...\n");
		}
		solve(trace_flag);
	}
	
	public void solve(boolean trace_flag) {
		simplex.dict.set_dual_view();
		while(true) {
			if(add_cutting_planes()) {
				iter++;
				try {
					simplex.solve(false);
				} catch (Exception e) {
					throw new ArithmeticException("ILP is INFEASIBLE!");
				}
				if(trace_flag) {
					StdOut.println(String.format(java.util.Locale.UK, "\niter %3d: current objective value: %.4f", iter, simplex.dict.z[0]));
				}
			} else {
				if(trace_flag) {
					print_solution();
				}
				break;
			}
		}
	}
	
	public boolean add_cutting_planes() {
		boolean is_added = false;
		int m = simplex.dict.m; //use old m value (during addition of cutting planes m increases);
		int n = simplex.dict.n;
		for(int i = 0; i < m; i++) {
			if(!is_integral(simplex.dict.b(i))) {
				is_added = true;
				double[] a = new double[simplex.dict.n];
				double b = -frac(simplex.dict.b(i));
				for(int j = 0; j < n; j++) {
					a[j] = frac(-simplex.dict.A(i, j));
				}
				simplex.dict.add_cutting_plane(a, b);
			}
		}
		return is_added;
	}
	
	public double[] get_solution() {
		return simplex.solution;
	}
	
	public void print_solution() {
		StdOut.println("\n*** Results ***\n");
		StdOut.println("Optimal solution obtained after " + iter + " cutting plane iterations" + ":\n");
		for(int i = 0; i < simplex.dict.n; i++) {
			StdOut.println(String.format(java.util.Locale.UK, "var " + (i+1) + ": %.4f", simplex.solution[i]));
		}
		StdOut.println(String.format(java.util.Locale.UK, "\nOptimal objective value: %.4f", simplex.dict.z[0]));
	}	
	
	public static double frac(double a) {
		return a - Math.floor(a);
	}
	
	public static boolean is_integral(double x) {
		return (frac(x) < tol) || (frac(x) > (1 - tol));
	}
	
	public static void main(String[] args) throws IOException {
		Stopwatch sw = new Stopwatch();
		boolean trace_flag = false;
		if("-t".equals(args[1])) {
			trace_flag = Boolean.parseBoolean(args[2]);
		}
		ILPSolver solver = new ILPSolver(new Dictionary(args[0]), trace_flag);
		StdOut.println("\nTiming results: " + sw.elapsedTime());
	}
}