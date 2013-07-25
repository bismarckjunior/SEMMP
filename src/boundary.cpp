#include "boundary.h" 

void setBoundaryConditions(int *j, int *Map, double *Ax, double *b, 
									 double N, double E, double W, double S, 
									 Boundary* boundary, Block *block)
{
	int bool_E, bool_W, bool_N, bool_S, i, k, m, tmp;
	int diag_index = *j;

	bool_E = bool_W = bool_N = bool_S = TRUE;
	
	setprogname("boundary conditions");

	for(k=block->nBoundary-1; k >= 0; k--){

		i = block->boundary[k];
		
		if (boundary[i].type == PRESSURE_ESPECIFIED){
			if(boundary[i].side == EAST && bool_E == TRUE){
				Ax[Map[diag_index]] -= E;
				*b -= 2*E*boundary[i].value;
				bool_E = FALSE;
			}
			else if(boundary[i].side == WEST && bool_W == TRUE){
				Ax[Map[diag_index]] -= W;
				*b -= 2*W*boundary[i].value;
				bool_W = FALSE;
			}
			else if(boundary[i].side == NORTH && bool_N == TRUE){
				Ax[Map[diag_index]] -= N;
				*b -= 2*N*boundary[i].value;
				bool_N = FALSE;
			}
			else if(boundary[i].side == SOUTH && bool_S == TRUE){
				Ax[Map[diag_index]] -= S;
				*b -= 2*S*boundary[i].value;
				bool_S = FALSE;
			}
			else{
				weprintf(
"there are more than one boundary conditions to the same\
                         block. Forgeting [Boundary %i] for block (%d,%d).", 
					i, block->row, block->col);	
				for(m=i; m < block->nBoundary-1; m++){
					tmp = block->boundary[m+1];
					block->boundary[m+1] = block->boundary[m];
					block->boundary[m] = tmp;
				}
				block->nBoundary--;
			}
		}
	
		else if (boundary[i].type == PRESSURE_GRADIENT_ESPECIFIED){
			if(boundary[i].side == EAST && bool_E == TRUE){
				Ax[Map[diag_index]] += E;
				*b += block->dx*E*boundary[i].value;
				bool_E = FALSE;
			}
			else if(boundary[i].side == WEST && bool_W == TRUE){
				Ax[Map[diag_index]] += W;
				*b -= block->dx*W*boundary[i].value;
				bool_W = FALSE;
			}
			else if(boundary[i].side == NORTH && bool_N == TRUE){
				Ax[Map[diag_index]] += N;
				*b -= block->dy*N*boundary[i].value;
				bool_N = FALSE;
			}
			else if(boundary[i].side == SOUTH && bool_S == TRUE){
				Ax[Map[diag_index]] += S;
				*b -= block->dy*S*boundary[i].value;
				bool_S = FALSE;
			}
			else{
				weprintf(
"there are more than one boundary conditions to the same\
                         block. Forgeting [Boundary %i] for block (%d,%d).", 
				  i, block->row, block->col);	
				for(m=i; m < block->nBoundary-1; m++){
					tmp = block->boundary[m+1];
					block->boundary[m+1] = block->boundary[m];
					block->boundary[m] = tmp;
				}
				block->nBoundary--;
			}
		}
	}

	(*j)++;
	if(E > 0.0 && bool_E == TRUE){
		Ax[Map[*j]] = E;
		(*j)++;
	}
	if(W > 0.0 && bool_W == TRUE){
		Ax[Map[*j]] = W;
		(*j)++;
	}
	if(N > 0.0 && bool_N == TRUE){
		Ax[Map[*j]] = N;
		(*j)++;
	}
	if(S > 0.0 && bool_S == TRUE){
		Ax[Map[*j]] = S;
		(*j)++;
	}

	unsetprogname();
}