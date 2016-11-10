Information for #404

The conceptual model presented to the user is described in 
	pandeia/ui/doc/ui_conceptual_model.txt
This document appears to be missing a bit here and there, but I don't
know what.

The client->server API is described in
	pandeia/ui/doc/transactions/ws_trans.txt
	pandeia/ui/doc/transactions/plot_trans.*

The client->engine API is not documented.
	-> try to find this

The server->engine API is documented in 
	pandeia/engine/doc/engine_api.rst
	pandeia/engine/doc/engine_input_api.rst
	pandeia/engine/doc/engine_output_api.rst
	pandeia/engine/doc/oapi_coords.fig
	pandeia/engine/doc/oapi_coords.gif

The server->database API is the Pandokia portable DBAPI wrapper.

The database schema is:
	pandeia/ui/pandeia/server/database/schema.sql
		(data relating to users)
	pandeia/ui/pandeia/server/starter/starter.sql
		(data relating to system startup)

	pandeia/ui/pandeia/server/database/workbooks.pdf
		(relation diagram of user data)

Initial workbooks are:
	pandeia/ui/pandeia/server/database/*.wb

client->database->engine API presentation in database.pptx , database.pdf
