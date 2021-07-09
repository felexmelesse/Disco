
#include "../paul.h"
#include "../boundary.h"

void boundary_trans( struct domain * theDomain , int dim )
{
    int nq = 1;
    int q[1];
    q[0] = UPP;

    if(dim == 1)
    {
        boundary_zerograd_rinn(theDomain, 1);
        boundary_fixed_q_rinn(theDomain, q, nq);
        boundary_fixed_rout(theDomain);
    }
    else if(dim == 2)
    {
        boundary_zerograd_zbot(theDomain, 1);
        boundary_zerograd_ztop(theDomain, 1);
    }
}

