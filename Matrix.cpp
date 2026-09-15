#include <stdio.h>
#include <stdlib.h>  
#include <math.h> 

double mu2(double x){
	return x * x; 
} 

void data_to_array(double **a, int m, int n, FILE *f){
	if(f == NULL){
		printf("Khong the mo file !\n"); 
	} 
	else{
		int i = 0, j = 0; 
		double p; 
		while(fscanf(f, "%lf", &p) != -1){
			for(int k = 0; i < m; k++){
				if(i % m == m - 1){
					a[m - 1][j] = p;
					j++; 
				} 
				else if(i % m == k){
					a[k][j] = p;
				} 
			} 
			i++;  
		} 
	} 
} 

double l2_norm(double *x, int n){
	double sum = 0;
	for(int i = 0; i < n; i++){
		sum += mu2(x[i]); 
	} 
	return sqrt(sum); 
} 

double distance(double *a, double *b, int n){
	double sum = 0;
	for(int i = 0; i < n; i++){
		sum += mu2(a[i] - b[i]); 
	} 
	return sqrt(sum); 
} 

int sign(double x){
	int s = x >= 0 ? 1 : -1; 
	return s; 
} 

double *scale(double *x, int n){
	double *x_scale = (double*)malloc(n * sizeof(double)); 
	double l2norm = l2_norm(x, n);
	if(l2norm == 0){
		return x; 
	} 
	for(int i = 0; i < n; i++){
		x_scale[i] = x[i] / l2norm;
	} 
	return x_scale; 
} 

void enter(double **a, int m, int n){ 
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			scanf("%lf", &a[i][j]); 
		}  
	}
} 

void print(double **a, int m, int n){
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			printf("%10lf ", a[i][j]); 
		} 
		printf("\n"); 
	} 
} 

double **matrix(int m, int n){
	double **a = (double**)malloc(m * sizeof(double*));
	for(int i = 0; i < m; i++){
		a[i] = (double*)malloc(n * sizeof(double)); 
	}
	return a; 
}

double **square_matrix(int n){
	double **a = (double**)malloc(n * sizeof(double*));
	for(int i = 0; i < n; i++){
		a[i] = (double*)malloc(n * sizeof(double)); 
	}
	return a; 
}

double **identity_matrix(int n){
	double **I = square_matrix(n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(i == j){
				I[i][j] = 1; 
			} 
			else{
				I[i][j] = 0;
			} 
		} 
	} 
	return I; 
} 


double **zero_matrix(int m, int n){
	double **Z = matrix(m, n);
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			Z[i][j] = 0; 
		} 
	} 
	return Z; 
} 

double **one_col_matrix(int n){
	double **one_col = matrix(n, 1);
	for(int i = 0; i < n; i++){
		one_col[i][0] = 1; 
	}  
	return one_col; 
} 

double **transpose(double **a, int m, int n){
	double **aT = matrix(n, m);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			aT[i][j] = a[j][i]; 
		} 
	} 
	return aT; 
} 

double **matrix_add(double **a, int m, int n, double **b){
	double **c = matrix(m, n);
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			c[i][j] = a[i][j] + b[i][j]; 
		} 
	} 
	return c; 
} 

double **ne_matrix(double **a, int m, int n){
	double **X = matrix(m, n); 
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			X[i][j] = -a[i][j]; 
		} 
	} 
	return X; 
} 

int check_equal(double **a, int m, int n, double **b){
	int cnt = 0;
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			if(a[i][j] != b[i][j]){
				cnt++; 
			} 
		} 
	} 
	if(cnt != 0){
		return 0; 
	} 
	return 1; 
} 

double **scalar_multi(double **a, int m, int n, double x){
	double **X = matrix(m, n);  
	for(int i = 0; i < m; i++){
		for(int j = 0; j < m; j++){
			X[i][j] = x * a[i][j]; 
		} 
	} 
	return X; 
} 

double **matrix_multi(double **a, int m, int n, double **b, int q){
	double **c = matrix(m, q);
	c = zero_matrix(m, q); 
	for(int i = 0; i < m; i++){
		for(int j = 0; j < q; j++){
			for(int k = 0; k < n; k++){
				c[i][j] += a[i][k] * b[k][j]; 
			} 
		} 
	} 
	return c; 
} 

double dot_product(double **a, double **b, int m, int n){
	double c = 0; 
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			c += a[i][j] * b[i][j];  
		} 
	} 
	return c; 
} 

double **Schur_product(double **a, double **b, int m, int n){
	double **c = matrix(m, n);
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			c[i][j] = a[i][j] * b[i][j];  
		} 
	} 
	return c;  
} 

double **matrix_pq(double **a, int m, int n, int p, int q){
	double **X = matrix(m - 1, n - 1);
	for(int i = 0; i < p - 1; i++){                 
		for(int j = 0; j < q - 1; j++){             
			X[i][j] = a[i][j];                     
		}                                           
	}                                                
	for(int i = 0; i < p - 1; i++){                  
		for(int j = q - 1; j < n - 1; j++){         
			X[i][j] = a[i][j + 1];                  
		}                                           
	}                                                
	for(int i = p - 1; i < m - 1; i++){              
		for(int j = 0; j < q - 1; j++){             
			X[i][j] = a[i + 1][j];                   
		}                                         
	} 
	for(int i = p - 1; i < m - 1; i++){
		for(int j = q - 1; j < n - 1; j++){
			X[i][j] = a[i + 1][j + 1]; 
		} 
	}  
	return X; 
} 

double **scale_matrix(double **a, int m, int n){
	double **multi = matrix_multi(a, m, n, one_col_matrix(n), 1); 
	double **mean_a = scalar_multi(multi, m, 1, 1.0 / n); 
	double **one_row_matrix = transpose(one_col_matrix(n), n, 1); 
	double **mean_d = matrix_multi(mean_a, m, 1, one_row_matrix, n);  
	double **scale_a = matrix_add(a, m, n, ne_matrix(mean_d, m, n)); 
	return scale_a; 
} 

double det(double **, int);

double cofactor(double **a, int n, int p, int q){
	int c = ((p + q) % 2 == 0) ? 1 : -1; 
	double Apq = c * det(matrix_pq(a, n, n, p, q), n - 1);
	return Apq; 
} 

double det(double **a, int n){
	if(n == 1){
		return a[0][0]; 
	} 
	double det = 0;
	for(int i = 0; i < n; i++){
		det += a[0][i] * cofactor(a, n, 1, i + 1); 
	} 
	return det; 
} 

double trace(double **a, int n){
	double x = 0;
	for(int i = 0; i < n; i++){
		x += a[i][i]; 
	} 
	return x; 
} 

double **inverse_matrix(double **a, int n){
	double det_a = det(a, n); 
	if(det_a == 0){
		printf("Khong ton tai ma tran nghich dao do dinh thuc bang 0\n");
		return zero_matrix(n, n); 
	} 
	double **X = square_matrix(n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			X[i][j] = cofactor(a, n, i + 1, j + 1); 
		} 
	} 
	X = scalar_multi(transpose(X, n, n), n, n, 1 / det_a);
	return X; 
} 

double Frobenius(double **a, int m, int n){
	double sum = 0;
	for(int i = 0; i < m; i++){
		for(int j = 0; j < n; j++){
			sum += mu2(a[i][j]); 
		} 
	} 
	return sum; 
}

double **Covariance(double **a, int m, int n){
	double **covariance = matrix(m, m);
	double **aT = transpose(a, m, n);
	covariance = scalar_multi(matrix_multi(a, m, n, aT, m), m, m, 1.0 / n); 
	return covariance; 
} 

double **Hessenberg(double **a, int n){
	for(int k = 1; k <= n - 2; k++){
		double *x = (double*)malloc((n - k) * sizeof(double));
		for(int i = 0; i < n - k; i++){
			x[i] = a[k + i][k - 1]; 
		} 
		double l2norm = l2_norm(x, n - k);
		double *v = (double*)malloc((n - k) * sizeof(double));
		double *e = (double*)malloc((n - k) * sizeof(double)); 
		for(int i = 0; i < n - k; i++){
			if(i == 0){
				e[i] = 1; 
			} 
			else{ 
			    e[i] = 0; 
			} 
		} 
		for(int i = 0; i < n - k; i++){
			v[i] = x[i] + sign(x[0]) * l2norm * e[i]; 
		} 
		v = scale(v, n - k);
		double **u = matrix(n - k, 1);
		for(int i = 0; i < n - k; i++){
			u[i][0] = v[i]; 
		} 
		double **P_k = matrix(n - k, n - k);
		double **uT = transpose(u, n - k, 1);
		double **uuT = matrix_multi(u, n - k, 1, uT, n - k);
		double **u2uT = scalar_multi(uuT, n - k, n - k, 2);
		P_k = matrix_add(identity_matrix(n - k), n - k, n - k, ne_matrix(u2uT, n - k, n - k));
		double **Pk = square_matrix(n);
		for(int i = 0; i < k; i++){
			for(int j = 0; j < k; j++){
				if(i == j){
					Pk[i][j] = 1; 
				} 
				else{
					Pk[i][j] = 0; 
				} 
			} 
		} 
		for(int i = 0; i < k; i++){
			for(int j = k; j < n; j++){
				Pk[i][j] = 0;  
			} 
		}
		for(int i = k; i < n; i++){
			for(int j = 0; j < k; j++){
				Pk[i][j] = 0;  
			} 
		}
		for(int i = k; i < n; i++){
			for(int j = k; j < n; j++){
				Pk[i][j] = P_k[i - k][j - k];  
			} 
		}
		a = matrix_multi(matrix_multi(Pk, n, n, a, n), n, n, transpose(Pk, n, n), n);  
	} 
	return a; 
}  

double uK(double **u){
	double delta = mu2(u[0][0] + u[1][1]) - 4 * (u[0][0] * u[1][1] - u[0][1] * u[1][0]);
	if(delta < 0){
		printf("Khong co tri rieng thuc, chi co tri rieng phuc !\n");
		return -3618; 
	} 
	double lamda1 = (u[0][0] + u[1][1] + sqrt(delta)) / 2; 
	double lamda2 = (u[0][0] + u[1][1] - sqrt(delta)) / 2;
	double uk = fabs(lamda1 - u[1][1]) < fabs(lamda2 - u[1][1]) ? lamda1 : lamda2;
	return uk; 
} 

double **Givens_matrix(int n, int i, int j, double c, double s){
	double **givens = identity_matrix(n); 
	givens[i][i] = c;
	givens[j][j] = c;
	givens[i][j] = s;
	givens[j][i] = -s; 
	return givens; 
} 

double ***Givens(double **a, int n){
	double **R = square_matrix(n);
	for(int i = 0; i < n; i ++){
		for(int j = 0; j < n; j++){
			R[i][j] = a[i][j]; 
		} 
	} 
	double **Q = identity_matrix(n);
	for(int i = 0; i < n - 1; i++){
		double c, s; 
		if(a[i + 1][i] == 0){
			c = 1, s = 0; 
		} 
		else if(fabs(a[i + 1][i]) > fabs(a[i][i])){
			double t = a[i][i] / a[i + 1][i];
			s = 1.0 / sqrt(1 + t * t);
			c = s * t; 
		} 
		else if(fabs(a[i][i]) >= fabs(a[i + 1][i])){
			double t = a[i + 1][i] / a[i][i];
			c = 1.0 / sqrt(1 + t * t);
			s = c * t; 
		} 
		double **givens = Givens_matrix(n, i, i + 1, c, s); 
		R = matrix_multi(givens, n, n, R, n);
		Q = matrix_multi(givens, n, n, Q, n); 
	}  
	Q = transpose(Q, n, n); 
	double ***QR = (double***)malloc(2 * sizeof(double**));
	QR[0] = Q;
	QR[1] = R; 
	return QR; 
} 

double ***LU(double **a, int n){
	double **L = identity_matrix(n);
	double **U = zero_matrix(n, n);
	for(int k = 0; k < n; k++){
		for(int j = k; j < n; j++){
			double sum = 0; 
			for(int m = 0; m <= k - 1; m++){
				sum += L[k][m] * U[m][j]; 
			} 
			U[k][j] = a[k][j] - sum; 
		}
		for(int i = k + 1; i < n; i++){
			double sum = 0;
			for(int m = 0; m <= k - 1; m++){
				sum += L[i][m] * U[m][k]; 
			} 
			L[i][k] = (a[i][k] - sum) / U[k][k];  
		} 
	} 
	double ***LU = (double***)malloc(2 * sizeof(double**));
	LU[0] = L;
	LU[1] = U;
	return LU; 
} 

double *AX0(double **a, int n, double *b){
	double ***LU_A = LU(a, n);
	double **L = LU_A[0];
	double **U = LU_A[1];
	double *x = (double*)malloc(n * sizeof(double)); 
	double *y = (double*)malloc(n * sizeof(double));
	y[0] = b[0]; 
	for(int i = 1; i < n; i++){ 
		double sum = 0; 
		for(int j = 0; j < i; j++){
			sum += L[i][j] * y[j]; 
		} 
		y[i] = b[i] - sum; 
	} 
	x[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for(int i = n - 2; i >= 0; i--){ 
		double sum = 0; 
		for(int j = i + 1; j < n; j++){
			sum += U[i][j] * x[j]; 
		} 
		x[i] = (y[i] - sum) / U[i][i]; 
	}
	return x; 
} 

double *Eigenvalue(double **a, int n){
	double *tri_r = (double*)malloc(n * sizeof(double)); 
	for(int i = 0; i < n; i++){
		tri_r[i] = 0; 
	} 
	while(1){
		double **u = square_matrix(2);
		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){
				u[i][j] = a[i + n - 2][j + n - 2];
			} 
		} 
		double uk = uK(u); 
		if(uk == -3618){
			printf("Khong the tinh toan voi tri rieng phuc!\n"); 
			return tri_r; 
		} 
		double **uKI = scalar_multi(identity_matrix(n), n, n, uk); 
		double **A = matrix_add(a, n, n, ne_matrix(uKI, n, n)); 
		double **H = Hessenberg(A, n);
		double ***QR = Givens(H, n); 
		double **RQ = matrix_multi(QR[1], n, n, QR[0], n);              
		a = matrix_add(RQ, n, n, uKI);  
		if(a[n - 1][n - 2] < 0.001){                                
			break; 
		} 
	}
	for(int i = 0; i < n; i++){
		tri_r[i] = a[i][i]; 
	} 
	return tri_r; 
} 

double rank_square_matrix(double **a, int n){
	double *x = Eigenvalue(a, n);
	double rank = 0; 
	for(int i = 0; i < n; i++){
		if(x[i] != 0){
			rank++; 
		} 
	} 
	return rank; 
} 

double **LamdaI(int n, int lamda){
	double **I = identity_matrix(n);
	return scalar_multi(I, n, n, lamda); 
}

double *eigenvector(double **a, int n, double r){
	double **A = matrix_add(a, n, n, ne_matrix(LamdaI(n, r), n, n));  
	double *xo = (double*)malloc(n * sizeof(double));
    for(int i = 1; i < n; i++){
		xo[i] = 1.0 / sqrt(n); 
	} 
	double *y = (double*)malloc(n * sizeof(double));  
	while(1){
		y = AX0(A, n, xo); 
		y = scale(y, n); 
		if(distance(xo, y, n) < 1){
			break; 
		} 
		xo = y;  
	} 
	return y; 
}	 

double *add_e(double *a, int n, double e){
	double *c = (double*)malloc(n * sizeof(double)); 
	for(int i = 0; i < n; i++){
		c[i] = a[i] + e; 
	} 
	return c; 
} 

double **Eigenvector(double **a, int n){
	double **vector = square_matrix(n); 
	double *tri_rieng = Eigenvalue(a, n);
	double *tri_r_e = add_e(tri_rieng, n, 0.0001); 
	for(int i = 0; i < n; i++){
		vector[i] = eigenvector(a, n, tri_r_e[i]); 
	}  
	return vector; 
} 

int main(){
	double **data_2d = matrix(2, 500);
	double **data_3d = matrix(3, 500);
	double **cov_2d = Covariance(data_2d, 2, 500);    //2x2 
	double **cov_3d = Covariance(data_3d, 3, 500);    //3x3 
//	int K;
//	scanf("%d", &K); 
	double **a = square_matrix(3);
	enter(a, 3, 3);
	 
//	print(inverse_matrix(a, 3), 3, 3); 
	
//	free(a); 
//	printf("%lf", det(a, 3));  
//	double *b = (double*)malloc(3 * sizeof(double));
//	for(int i = 0; i < 3; i++){
//		scanf("%lf", &b[i]); 
//	} 
	double *x = Eigenvalue(a, 3);
	for(int i = 0; i < 3; i++){
		printf("%lf ", x[i]); 
	}  
//	double *vtr = eigenvector(a, 3, 1.01); 
//	double **vtr = Eigenvector(a, 3);
} 
