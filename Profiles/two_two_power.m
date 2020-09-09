function data_out = two_two_power(s, am_in)
%
a0 = am_in(1);
a1 = am_in(2);
a2 = am_in(3);
a3 = am_in(4);
a4 = am_in(5);
a5 = am_in(6);

data_out = a0 * (  (   a3   * (1 - s.^a1).^a2)  +  ...
    ( (1-a3) * (1 - s.^a4).^a5)    ) ;

