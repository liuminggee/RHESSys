/*--------------------------------------------------------------*/
/*                                                              */
/*        sort_patch_layers                                    */
/*                                                              */
/*  NAME                                                        */
/*        sort_patch_layers                                    */
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  sort_patch_layers( struct patch_object *patch)                    */
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*    sorts canopy_stratum within a patch by height into    */
/*    different layers                    */
/*                                */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include "rhessys.h"
#include "functions.h"

void sort_patch_layers( struct patch_object *patch)
{
    /*--------------------------------------------------------------*/
    /*  Local function declaration                                  */
    /*--------------------------------------------------------------*/
    int key_compare(const void * e1,const  void *e2 );
    void    *alloc(     size_t, char *, char *);
    /*--------------------------------------------------------------*/
    /*  Local variable definition.                                  */
    /*--------------------------------------------------------------*/
    int i, j,k;
    int list_bottom;
    double cover_fraction;
    /*--------------------------------------------------------------*/
    /*    free current layer structure                */
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*        Establish index of next free list entry.    */
    /*--------------------------------------------------------------*/

    //10062025LML
    //1# solve the problem for total canopy cover > 1 for one layer
    //2# making sure the tree strata is higher than understory strata
    //Note: it's a simple fix!!!

    bool need_height_adj_for_total_cover_fraction = false;
    bool need_height_adj_for_overstory_tree = false;

    do { //total cover fraction of each layer cannot over 1.0
        for ( i=0 ; i<patch[0].num_canopy_strata; i++ ) {
            free(patch[0].layers[i].strata);
        }
        do { //overstory tree need to be on the top, otherwise it may not grow under the shade. 10062025 LML
            need_height_adj_for_overstory_tree = false;
            for ( i=0 ; i<patch[0].num_canopy_strata; i++ ) {
                patch[0].layers[i].height = 0;  //10062025LML
                patch[0].layers[i].count = 0;
            }
            patch[0].num_layers = 0;
            list_bottom = 0;
            /*--------------------------------------------------------------*/
            /*    Determine the unique height layers in the patch    */
            /*--------------------------------------------------------------*/
            for( i=0; i<patch[0].num_canopy_strata ; i++ ){
                /*--------------------------------------------------------------*/
                /*        Check if this height alread exists.        */
                /*--------------------------------------------------------------*/
                j = 0;
                while( (j<list_bottom) && ! close_enough(patch[0].canopy_strata[i][0].epv.height,
                    patch[0].layers[j].height) ){
                    j++;
                }
                /*--------------------------------------------------------------*/
                /*        If we did not find this height in the list    */
                /*        add it to the bottom of the list.        */
                /*        otherwise just increment the count of layers    */
                /*        at height "height".                */
                /*--------------------------------------------------------------*/
                if ( j >= list_bottom ){
                    (patch[0].layers[list_bottom]).height =
                        patch[0].canopy_strata[i][0].epv.height;
                    (patch[0].layers[list_bottom]).count = 1;
                    list_bottom++;
                }
                else {
                    (patch[0].layers[j]).count++;
                    //10102023LML note: if patch and strata height are zero, multiple
                    //strates will be in this zero layer, i.e. count > 1 while the total
                    //cover fraction might be higher than 1
                }
            }
            /*--------------------------------------------------------------*/
            /*    Define number of unique layers in this patch.        */
            /*--------------------------------------------------------------*/
            patch[0].num_layers = list_bottom;
            /*--------------------------------------------------------------*/
            /*    Now sort the layer list into descending order.        */
            /*--------------------------------------------------------------*/
            qsort(
                (void *) patch[0].layers,
                (size_t) patch[0].num_layers,
                sizeof(struct layer_object),
                &key_compare);
            /*--------------------------------------------------------------*/
            /*    Now construct a list of pointers to strata at each    */
            /*    height layer                        */
            /*--------------------------------------------------------------*/
            //10062025LML make sure tree strata is the highest
            double top_layer_height = patch[0].layers[0].height;
            for ( j=0 ; j<patch[0].num_canopy_strata; j++ ){
                if ((patch[0].canopy_strata[j][0].defaults[0][0].epc.veg_type == TREE)
                    && (patch[0].canopy_strata[j][0].defaults[0][0].ID < 10) //overstory strata
                    && ((patch[0].canopy_strata[j][0].epv.height < top_layer_height)
                        || close_enough(top_layer_height,0.0))) {
                    patch[0].canopy_strata[j][0].epv.height = top_layer_height + 0.001;
                    need_height_adj_for_overstory_tree = true;
                }
            }
        } while (need_height_adj_for_overstory_tree);

        for ( i=0 ; i<patch[0].num_layers ; i++ ){
            /*--------------------------------------------------------------*/
            /*        Allocate the list for layer i            */
            /*--------------------------------------------------------------*/
            patch[0].layers[i].strata = (long *)
                alloc(patch[0].layers[i].count*sizeof(long),
                "patch[0].layers[i].strata",
                "construct_patch");
            /*--------------------------------------------------------------*/
            /*        Reset the cover_fraction accumulator        */
            /*--------------------------------------------------------------*/
            cover_fraction = 0.0;
            /*--------------------------------------------------------------*/
            /*    assign a bottom of layer                */
            /*--------------------------------------------------------------*/
            if (i != patch[0].num_layers - 1)
                patch[0].layers[i].base = patch[0].layers[i+1].height;
            else
                patch[0].layers[i].base = 0.0;
            /*--------------------------------------------------------------*/
            /*        Find all strata with height matching layer i    */
            /*--------------------------------------------------------------*/
            //10102023LML
            double max_cover_fraction = 0;
            for ( j=0 ; j<patch[0].num_canopy_strata; j++ ){
                k = 0;
                /*--------------------------------------------------------------*/
                /*            check if this stratum has layer i height*/
                /*--------------------------------------------------------------*/
                //06122023LML note: if layers height is zero and epv height is zero
                //the total cover might greater than one, which will affect strata process and mass balance!
                if (close_enough(patch[0].canopy_strata[j][0].epv.height,
                    patch[0].layers[i].height)) {
                /*--------------------------------------------------------------*/
                /*        Add the stratum index to the layer if it matches*/
                /*--------------------------------------------------------------*/
                    patch[0].layers[i].strata[k] = j;
                    /*--------------------------------------------------------------*/
                    /*        Keep a running total of the cover fraction in    */
                    /*        this layer to check that it adds to 1.0        */
                    /*--------------------------------------------------------------*/
                    if (patch[0].canopy_strata[j][0].cover_fraction > max_cover_fraction)
                        max_cover_fraction = patch[0].canopy_strata[j][0].cover_fraction; //10102023LML
                    cover_fraction += patch[0].canopy_strata[j][0].cover_fraction;
                    k++;
                }
            }
        /*--------------------------------------------------------------*/
        /*        Report a fatal error if the cover fraction for    */
        /*        this layer does not add to 1.0            */
        /*--------------------------------------------------------------*/
            if ( cover_fraction > 1.0 ){
                printf( "\nWARNING: in sort_patch_layers cover fraction of layer height %f greater than 1.0! \nAdjustment...\n"
                        ,cover_fraction);
                need_height_adj_for_overstory_tree = true;
                //increase the height for trees
                for (j = 0; j < patch[0].layers[i].count; j++) {
                    int strata_idx = patch[0].layers[i].strata[j];
                    if (patch[0].canopy_strata[strata_idx][0].defaults[0][0].epc.veg_type == TREE) {
                        patch[0].canopy_strata[strata_idx][0].epv.height += 0.0001 * (1 + j);
                    }
                }
            }
            patch[0].layers[i].null_cover = 1.0 - max_cover_fraction;
        } //layer i
    } while (need_height_adj_for_overstory_tree);
    return;
}
