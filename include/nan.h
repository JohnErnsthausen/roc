/* 
© 2002-2017.  Nedialko S. Nedialkov, McMaster University, John
D. Pryce, Cardiff University.  All rights reserved.

The DAETS package may be downloaded and used for personal, research,
academic and non-commercial purposes only, and all other uses are
strictly prohibited.  If you wish to download and use the DAETS
package for any other purpose, please contact nedialk@mcmaster.ca.
The DAETS package is provided “as is” without warranty of any kind,
expressed or implied.  In no event shall the authors or McMaster
University or Cardiff University be liable for any claim, damages or
liability, arising from, out of or in connection with DAETS or the use
or other dealings in DAETS.
*/

#ifndef NAN_H
#define NAN_H

void init_with_nans(int n, double *d);
int check_for_nan(int n, double *x);

#endif // NAN_H
