adios_groupsize = 0;
adios_totalsize = 0;
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_buf_size = 4;
adios_read (adios_handle, "NX", &NX, adios_buf_size);
adios_buf_size = 8 * (NX);
adios_read (adios_handle, "/temperature", t, adios_buf_size);
