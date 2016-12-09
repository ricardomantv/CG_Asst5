// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
      const T& position0,
      const T& position1,
      const T& tangent0,
      const T& tangent1,
      double normalizedTime,
      int derivative )
{
   // TODO IMPLEMENT ME (TASK 1A)
   double t = normalizedTime;
   double t2 = pow(t, 2);
   double t3 = pow(t, 3);

   double h00, h10, h01, h11;
   if(derivative == 0) {
     h00 = 2 * t3 - 3 * t2 + 1;
     h10 = t3 - 2 * t2 + t;
     h01 = -2 * t3 + 3 * t2;
     h11 = t3 - t2;
   } else if(derivative == 1) {
     h00 = 6 * t2 - 6 * t;
     h10 = 3 * t2 - 4 * t + 1;
     h01 = -6 * t2 + 6 * t;
     h11 = 3 * t2 - 2 * t;
   } else if(derivative == 2) {
     h00 = 12 * t - 6;
     h10 = 6 * t - 4;
     h01 = -12 * t + 6;
     h11 = 6 * t - 2;
   } else {
     // Invalid derivative
     return T();
   }

   return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

// Returns a state interpolated between the values directly before and after the given time.
template <class T>
inline T Spline<T>::evaluate( double time, int derivative )
{
   // TODO IMPLEMENT ME (TASK 1B)
   if (knots.size() < 1) return T();

   if(knots.size() == 1) {
     if(derivative == 0) {
       return knots.begin()->second;
     } else {
       return T();
     }
   } else {
     if(time <= knots.begin()->first) {
       // Return first knot if query is before it
       if(derivative == 0) {
         return knots.begin()->second;
       } else {
         return T();
       }

     }

     if(prev(knots.end())->first <= time) {
       // Return last knot if query is after it
       if(derivative == 0) {
         return prev(knots.end())->second;
       } else {
         return T();
       }
     }

     // KnotIter k1_it = knots.lower_bound(time);
     KnotIter k2_it = knots.upper_bound(time);
     KnotIter k1_it = prev(k2_it);

     double t1 = k1_it->first;
     double t2 = k2_it->first;
     T k1 = k1_it->second;
     T k2 = k2_it->second;

     double t0, t3;
     T k0, k3;
     if(k1_it == knots.begin()) {
       // No second knot to left of queried time, create virtual knot
       t0 = t1 - (t2 - t1);
       k0 = k1 - (k2 - k1);
     } else {
       t0 = prev(k1_it)->first;
       k0 = prev(k1_it)->second;
     }

     if(k2_it == prev(knots.end())) {
       // No second knot to right of queried time, create virtual knot
       t3 = t2 + (t2 - t1);
       k3 = k2 + (k2 - k1);
     } else {
       t3 = next(k2_it)->first;
       k3 = next(k2_it)->second;
     }

     double norm_t0 = (t0 - t1) / (t2 - t1); // will be negative
     double norm_t1 = 0;
     double norm_t2 = 1;
     double norm_t3 = (t3 - t1) / (t2 - t1); // will be greater than 1

     T m1 = (k2 - k0) / (norm_t2 - norm_t0);
     T m2 = (k3 - k1) / (norm_t3 - norm_t1);

     double normT = (time - t1) / (t2 - t1);

     return cubicSplineUnitInterval(k1, k2, m1, m2, normT, derivative);
   }

}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance )
{
   // Empty maps have no knots.
   if( knots.size() < 1 )
   {
      return false;
   }

   // Look up the first element > or = to time.
   typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
   typename std::map<double, T>::iterator t1_iter;
   t1_iter = t2_iter;
   t1_iter--;

   if( t2_iter == knots.end() )
   {
      t2_iter = t1_iter;
   }

   // Handle tolerance bounds,
   // because we are working with floating point numbers.
   double t1 = (*t1_iter).first;
   double t2 = (*t2_iter).first;

   double d1 = fabs(t1 - time);
   double d2 = fabs(t2 - time);


   if(d1 < tolerance && d1 < d2)
   {
      knots.erase(t1_iter);
      return true;
   }

   if(d2 < tolerance && d2 < d1)
   {
      knots.erase(t2_iter);
      return t2;
   }

   return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue( double time, T value )
{
   knots[ time ] = value;
}

template <class T>
inline T Spline<T>::operator()( double time )
{
   return evaluate( time );
}
