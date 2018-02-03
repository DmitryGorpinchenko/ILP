import java.io.*;
import java.util.ArrayList;

class Dictionary {

	public int n;
	public int m;
	public ArrayList<Integer> basic = new ArrayList<>();
	private int[] non_basic;
	private ArrayList<Double> b = new ArrayList<>();
	private ArrayList<double[]> A = new ArrayList<>();
	public double[] z;
	public boolean is_final;
	public boolean is_unbounded;
	public boolean is_dual;
	public static double tol = 1e-6;
	
	public Dictionary(String file) throws IOException {
		parse(file);
	}
	
	public void set_dual_view() {
		is_dual = true;
		is_unbounded = false;
		is_final = false;
	}
	
	public void set_primal_view() {
		is_dual = false;
		is_unbounded = false;
		is_final = false;
	}
	
	public int non_basic(int id) {
		return is_dual ? basic.get(id) : non_basic[id];
	}
	
	public int basic(int id) {
		return is_dual ? non_basic[id] : basic.get(id);
	}
	
	public double A(int i, int j) {
		return A.get(i)[j];
	}
	
	public double b(int i) {
		return b.get(i);
	}
	
	public int n() {
		return is_dual ? m : n;
	}
	
	public int m() {
		return is_dual ? n : m;
	}

	public void add_cutting_plane(double[] a, double b) {
		m++;
		basic.add(n+m);
		this.b.add(b);
		A.add(a);
		is_final = false;
		is_unbounded = false;
	}	
	
	public double[] primal_solution() {
		double[] solution = new double[n];
		for(int i = 0; i < m; i++) {
			if(basic.get(i) <= n) {
				solution[basic.get(i)-1] = b.get(i);
			}
		}
		return solution;
	}
	
	public int get_entering() {
		if(is_dual) {
			return dual_entering();
		}
		return primal_entering();
	}
	
	public int get_leaving(int enter_id) {
		if(is_dual) {
			return dual_leaving(enter_id);
		}
		return primal_leaving(enter_id);
	}
	
	public void pivot(int enter_id, int leave_id) {
		if(is_dual) {
			dual_pivot(enter_id, leave_id);
		} else {
			primal_pivot(enter_id, leave_id);
		}
	}
	
	public int primal_entering() {
		int var_num = Integer.MAX_VALUE;
		int id = -1;
		for(int i = 1; i <= n; i++) {
			if(/*z[i] > 0*/ z[i] > tol && non_basic[i-1] < var_num) {
				var_num = non_basic[i-1];
				id = i-1;
			}
		}
		if(var_num == Integer.MAX_VALUE) {
			is_final = true;
		}
		return id;
	}
	
	public int primal_leaving(int enter_id) {
		int var_num = Integer.MAX_VALUE;
		int id = -1;
		double curr = Double.POSITIVE_INFINITY;
		for(int i = 0; i < m; i++) {
			if(/* A.get(i)[enter_id] < 0*/ (A.get(i)[enter_id] < -tol)) {
				double ratio = -b.get(i)/A.get(i)[enter_id];
				if(curr > ratio) {
					curr = ratio;
					var_num = basic.get(i);
					id = i;
				} else if(curr == ratio && var_num > basic.get(i)) {
					var_num = basic.get(i);
					id = i;
				} 
			}
		}
		if(var_num == Integer.MAX_VALUE) {
			is_unbounded = true;
		}
		return id;
	}
	
	public void primal_pivot(int enter_id, int leave_id) {
		//leave and enter to the basis
		int temp = non_basic[enter_id];
		non_basic[enter_id] = basic.get(leave_id);
		basic.set(leave_id, temp);
		//update leaving row in the dictionary
		b.set(leave_id, -b.get(leave_id)/A.get(leave_id)[enter_id]);
		for(int i = 0; i < n; i++) {
			if(i != enter_id) {
				A.get(leave_id)[i] /= -A.get(leave_id)[enter_id];
			} 
		}
		A.get(leave_id)[enter_id] = 1/A.get(leave_id)[enter_id];
		//update other rows
		for(int i = 0; i < m; i++) {
			if(i != leave_id) {
				b.set(i, b.get(i) + b.get(leave_id)*A.get(i)[enter_id]);
				for(int j = 0; j < n; j++) {
					if(j != enter_id) {
						A.get(i)[j] += A.get(i)[enter_id]*A.get(leave_id)[j];
					}
				}
				A.get(i)[enter_id] = A.get(i)[enter_id]*A.get(leave_id)[enter_id];
			}
		}
		//update objective coefficients
		z[0] += z[enter_id+1]*b.get(leave_id);
		for(int i = 1; i <= n; i++) {
			if(i-1 != enter_id) {
				z[i] += z[enter_id+1]*A.get(leave_id)[i-1];
			}
		}
		z[enter_id+1] = z[enter_id+1]*A.get(leave_id)[enter_id];
	}	
	
	public int dual_entering() {
		int var_num = Integer.MAX_VALUE;
		int id = -1;
		for(int i = 0; i < m; i++) {
			if(/* b.get(i) < 0 */ b.get(i) < -tol && basic.get(i) < var_num) {
				var_num = basic.get(i);
				id = i;
			}
		}
		if(var_num == Integer.MAX_VALUE) {
			is_final = true;
		}
		return id;
	}
	
	public int dual_leaving(int enter_id) {
		int var_num = Integer.MAX_VALUE;
		int id = -1;
		double curr = Double.POSITIVE_INFINITY;
		for(int i = 0; i < n; i++) {
			if(/* A.get(enter_id)[i] > 0 */ (A.get(enter_id)[i] > tol)) {
				double ratio = -z[i+1]/A.get(enter_id)[i];
				if(curr > ratio) {
					curr = ratio;
					var_num = non_basic[i];
					id = i;
				} else if(curr == ratio && var_num > non_basic[i]) {
					var_num = non_basic[i];
					id = i;
				} 
			}
		}
		if(var_num == Integer.MAX_VALUE) {
			is_unbounded = true;
		}
		return id;
	}
	
	public void dual_pivot(int enter_id, int leave_id) {
		//leave and enter to the basis
		int temp = basic.get(enter_id);
		basic.set(enter_id, non_basic[leave_id]);
		non_basic[leave_id] = temp;
		//update leaving row in the dictionary
		z[leave_id+1] = z[leave_id+1]/A.get(enter_id)[leave_id];
		for(int i = 0; i < m; i++) {
			if(i != enter_id) {
				A.get(i)[leave_id] /= A.get(enter_id)[leave_id];
			} 
		}
		A.get(enter_id)[leave_id] = 1/A.get(enter_id)[leave_id];
		//update other rows
		for(int i = 0; i < n; i++) {
			if(i != leave_id) {
				z[i+1] += -z[leave_id+1]*A.get(enter_id)[i];
				for(int j = 0; j < m; j++) {
					if(j != enter_id) {
						A.get(j)[i] += -A.get(enter_id)[i]*A.get(j)[leave_id];
					}
				}
				A.get(enter_id)[i] = -A.get(enter_id)[i]*A.get(enter_id)[leave_id];
			}
		}
		//update objective coefficients
		z[0] += -b.get(enter_id)*z[leave_id+1];
		for(int i = 0; i < m; i++) {
			if(i != enter_id) {
				b.set(i, b.get(i) - b.get(enter_id)*A.get(i)[leave_id]);
			}
		}
		b.set(enter_id, -b.get(enter_id)*A.get(enter_id)[leave_id]);
	}
	
	public void parse(String file) throws IOException {
		BufferedReader input = new BufferedReader(new FileReader(file));
		String line = input.readLine();
		String[] tokens = line.split("\\s+");
		m = Integer.parseInt(tokens[0]);
		n = Integer.parseInt(tokens[1]);
		line = input.readLine();
		tokens = line.split("\\s+");	
		for(int i = 0; i < m; i++) {
			basic.add(Integer.parseInt(tokens[i]));
		}
		non_basic = new int[n];
		line = input.readLine();
		tokens = line.split("\\s+");	
		for(int i = 0; i < n; i++) {
			non_basic[i] = Integer.parseInt(tokens[i]);
		}
		line = input.readLine();
		tokens = line.split("\\s+");	
		for(int i = 0; i < m; i++) {
			b.add(Double.parseDouble(tokens[i]));
		}
		for(int i = 0; i < m; i++) {
			line = input.readLine();
			tokens = line.split("\\s+");
			A.add(new double[n]);
			for(int j = 0; j < n; j++) {
				A.get(i)[j] = Double.parseDouble(tokens[j]);
			}
		}
		z = new double[n+1];
		line = input.readLine();
		tokens = line.split("\\s+");
		for(int i = 0; i < n+1; i++) {
			z[i] = Double.parseDouble(tokens[i]);
		}
		input.close();
	}
	
	public String toString() {
		String s = "";
		s += (m + " " + n + "\n");
		for(int i = 0; i < m; i++) {
			s += basic.get(i) + " ";
		}
		s += "\n";
		for(int i = 0; i < n; i++) {
			s += non_basic[i] + " ";
		}
		s += "\n";
		for(int i = 0; i < m; i++) {
			s += String.format(java.util.Locale.UK, "%.2f ", b.get(i));
		}
		s += "\n";
		for(int i = 0; i < m; i++) {
			for(int j = 0; j < n; j++) {
				s += String.format(java.util.Locale.UK, "%.2f ", A.get(i)[j]);
			}
			s += "\n";
		}
		for(int i = 0; i < n+1; i++) {
			s += String.format(java.util.Locale.UK, "%.2f ", z[i]);
		}
		return s;
	}
}