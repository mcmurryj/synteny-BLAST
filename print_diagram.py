0.  Parse parms.  Get the table, as well as the # of ortholologs required to map.
1.Read the table to a dataframe.  Filter dataframe to rows with NA NA AN NA or similar.

Return a dict:
keys = the WP IDs of the proteins from the cluster in question. Were Na Na NA appears in the row.
values = (start, stop, pfam IDs)
2.  For each WP ID in keys:
    Count dimensions of dataframe with all rows where cluster_member_ID = wp_id
    if count > cutoff, assign a color from color_list

3.  Get min and max of startstop; subtract min from all values

4.
