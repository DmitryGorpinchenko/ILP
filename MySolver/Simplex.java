import java.io.*;

class Simplex {
	
	public Dictionary dict;
	public double[] solution;
	public int iter = 0;
	public int enter_var;
	public int leave_var;
	
	public Simplex(Dictionary dict, boolean trace_flag) {
		this.dict = dict;
		if(initialize(trace_flag)) {
			solve(trace_flag);
		} else {
			throw new ArithmeticException("Linear program is INFEASIBLE!");
		}
	}
	
	//returns true in case of obtaining a feasible initial dictionary, false otherwise
	public boolean initialize(boolean trace_flag) {
		boolean is_needed = false;
		for(int i = 0; i < dict.m; i++) {
			if(dict.b(i) < 0) {
				is_needed = true;
				break;
			}
		}
		if(!is_needed) {
			return true;
		}
		try {
			if(trace_flag) {
				StdOut.println("\nInitialization phase ... ");
			}
			//copy an initial objective in order to restore it later on
			double[] primal_z = new double[dict.n + 1];
			for(int i = 0; i <= dict.n; i++) {
				primal_z[i] = dict.z[i];
			}
			//change objective in order to proceed with initialization
			for(int i = 1; i <= dict.n; i++) {
				dict.z[i] = -1.0;
			}
			//set the dual view of dictionary
			this.dict.set_dual_view();
			solve(false);
			//set back the primal view of dictionary
			this.dict.set_primal_view();
			//restore an objective
			this.dict.z = new double[this.dict.n+1];
			for(int i  = 0; i < this.dict.n; i++) {
				int non_basic_id = this.dict.non_basic(i);
				if(non_basic_id <= this.dict.n) {
					this.dict.z[i+1] = primal_z[non_basic_id];
				}
			}
			for(int i = 0; i < this.dict.m; i++) {
				int id = this.dict.basic(i);
				if(id <= this.dict.n) {
					this.dict.z[0] += primal_z[id]*this.dict.b(i);
					for(int j = 1; j <= this.dict.n; j++) {
						this.dict.z[j] += primal_z[id]*this.dict.A(i, j-1);
					}
				}
			}				
		} catch (ArithmeticException e) {
			return false;
		} 
		return true;
	}
	
	public void solve(boolean trace_flag) {
		iter = 0;
		if(trace_flag) {
			StdOut.println("\nOptimization phase ... ");
		}
		while(true) {
			int enter_id = dict.get_entering();
			if(dict.is_final) {
				solution = dict.primal_solution();
				if(trace_flag) {
					print_solution();
				}
				break;
			}
			int leave_id = dict.get_leaving(enter_id);
			if(dict.is_unbounded) {
				throw new ArithmeticException("Linear program is UNBOUNDED!");
			}
			iter++;
			enter_var = dict.non_basic(enter_id); //cache entering and leaving vars before it would be changed by pivot() 
			leave_var = dict.basic(leave_id);	   			
			dict.pivot(enter_id, leave_id);
			if(trace_flag) {
				print_trace_info();
			}
		}
	}
	
	public void print_trace_info() {
		StdOut.println("\nIteration " + iter + ":\n");
		StdOut.println(" - Entering var: " + enter_var);  
		StdOut.println(" - Leaving var: " + leave_var);
		StdOut.println(" - Current objective value: " + dict.z[0]);
	}
	
	public void print_solution() {
		StdOut.println("\nOptimal solution obtained after " + iter + " simplex iterations" + ":\n");
		for(int i = 0; i < dict.n; i++) {
			StdOut.println(String.format(java.util.Locale.UK, "var " + (i+1) + ": %.4f", solution[i]));
		}
		StdOut.println(String.format(java.util.Locale.UK, "\nOptimal objective value: %.4f", dict.z[0]));
	}
	
	public double[] get_solution() {
		return solution;
	}
	
	public static void main(String[] args) throws IOException {
		Stopwatch sw = new Stopwatch();
		boolean trace_flag = false;
		if("-t".equals(args[1])) {
			trace_flag = Boolean.parseBoolean(args[2]);
		}
		Simplex solver = new Simplex(new Dictionary(args[0]), trace_flag);
		StdOut.println("\nTiming results: " + sw.elapsedTime());
	}
}