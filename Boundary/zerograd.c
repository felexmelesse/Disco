
#include "../paul.h"
#include <string.h>

void boundary_fixed_rinn( struct domain *theDomain);
void boundary_fixed_rout( struct domain *theDomain);
void boundary_fixed_zbot( struct domain *theDomain);
void boundary_fixed_ztop( struct domain *theDomain);
void boundary_zerograd_rinn( struct domain *theDomain, int diode);
void boundary_zerograd_rout( struct domain *theDomain, int diode);
void boundary_zerograd_zbot( struct domain *theDomain, int diode);
void boundary_zerograd_ztop( struct domain *theDomain, int diode);
void boundary_reflect_rinn( struct domain *theDomain);
void boundary_reflect_rout( struct domain *theDomain);
void boundary_reflect_zbot( struct domain *theDomain);
void boundary_reflect_ztop( struct domain *theDomain);
void boundary_fixed_horizon( struct domain *theDomain);

void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_zerograd_rinn(theDomain, 0);
        boundary_zerograd_rout(theDomain, 0);
    }
    else if(dim == 2)
    {
        boundary_zerograd_zbot(theDomain, 0);
        boundary_zerograd_ztop(theDomain, 0);
    }
}

