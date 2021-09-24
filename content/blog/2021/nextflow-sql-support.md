title=Introducing Nextflow support for SQL databases
date=2021-09-16
type=post
tags=nextflow,plugins,sql
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

The recent tweet introducing the [Nextflow support for SQL databases](https://twitter.com/PaoloDiTommaso/status/1433120149888974854) raised a lot of positive reaction. In this post, I want to describe more in detail how this extension works. 


Nextflow was designed with the idea to streamline the deployment of complex data pipelines in a scalable, portable and reproducible manner across different computing platforms. To make this all possible, it was decided the resulting pipeline and the runtime should be self-contained i.e. to not depend on separate services such as database servers. 

This makes the resulting pipelines easier to configure,  deploy, and allows for testing them using [CI services](https://en.wikipedia.org/wiki/Continuous_integration), which is a critical best practice for delivering high-quality and stable software.

Another important consequence is that Nextflow pipelines do not retain the pipeline state on separate storage. Said in a different way, the idea was - and still is - to promote stateless pipeline execution in which the computed results are only determined by the pipeline inputs and the code itself, which is consistent with the *functional* dataflow paradigm on which Nextflow is based. 

However, the ability to access SQL data sources can be very useful in data pipelines, for example, to ingest input metadata or to store task executions logs. 


### How does it work?

The support for SQL databases in Nextflow is implemented as an optional plugin component. This plugin provides two new operations into your Nextflow script: 

1) `fromQuery` performs a SQL query against the specified database and returns a Nextflow channel emitting them. This channel can be used in your pipeline as any other Nextflow channel to trigger the process execution with the corresponding values. 
2) `sqlInsert` takes the values emitted by a Nextflow channel and inserts them into a database table. 

The plugin supports out-of-the-box popular database servers such as MySQL, PostgreSQL and MariaDB. It should be noted that the technology is based on the Java JDBC database standard, therefore it could easily support any database technology implementing a driver for this standard interface. 

Disclaimer: This plugin is a preview technology. Some features, syntax and configuration settings can change in future releases.

### Let's get started!

The use of the SQL plugin requires the use of Nextflow 21.08.0-edge or later. If are using an older version, check [this page](https://www.nextflow.io/docs/latest/getstarted.html#stable-edge-releases) on how to update to the latest edge release.

To enable the use of the database plugin, add the following snippet in your pipeline configuration file. 

```
plugins {
  id 'nf-sqldb@0.1.0'
}
```

It is then required to specify the connection *coordinates* of the database service you want to connect to in your pipeline. This is done by adding a snippet similar to the following in your configuration file: 

```
sql {
    db {
        'my-db' {
              url = 'jdbc:mysql://localhost:3306/demo'
              user = 'my-user'
              password = 'my-password'
            }
    }
}
``` 

In the above example, replace `my-db` with a name of your choice (this name will be used in the script to reference the corresponding database connection coordinates). Also, provide a `url`, `user` and `password` matching your database server. 

Your script should then look like the following: 

```
nextflow.enable.dsl=2

process myProcess {
  input:
    tuple val(sample_id), path(sample_in) 
  output:
    tuple val(sample_id), path('sample.out') 

  """
  your_command --input $sample_id > sample.out
  """
}

workflow {

  query = 'select SAMPLE_ID, SAMPLE_FILE from SAMPLES'
  channel.sql.fromQuery(query, db: 'my-db') \
    | myProcess \
    | sqlInsert(table: 'RESULTS', db: 'my-db')

}
```

The above example shows how to perform a simple database query, pipe the results to a fictitious process named `myProcess` and finally store the process outputs into a database table named `RESULTS`.  

It is worth noting that Nextflow allows the use of any number of database instances in your pipeline, simply defining them in the configuration file using the syntax shown above. This could be useful to fetch database data from one data source and store the results into a different one. 

Also, this makes it straightforward to write [ETL](https://en.wikipedia.org/wiki/Extract,_transform,_load) scripts that span across multiple data sources.   

Find more details about the SQL plugin for Nextflow at [this link](https://github.com/nextflow-io/nf-sqldb).

## What about the self-contained property?  

You may wonder if adding this capability breaks the self-contained property of Nextflow pipelines which allows them to be run in a single command and to be tested with continuous integration services e.g. GitHub Action.   

The good news is that it does not ... or at least it should not if used properly. 

In fact, the SQL plugin includes the [H2](http://www.h2database.com/html/features.html) embedded in-memory SQL database that is used by default when no other database is provided in the Nextflow configuration file and can be used for developing and testing your pipeline without the need for a separate database service.  

Tip: Other than this, H2 also provides the capability to access and query CSV/TSV files as SQL tables. Read more about this feature at [this link](http://www.h2database.com/html/tutorial.html?highlight=csv&search=csv#csv).

### Conclusion

The use of this plugin adds to Nextflow the capability to query and store data into the SQL databases. Currently, the most popular SQL technologies are supported such as MySQL, PostgreSQL and MariaDB. In the future, support for other database technologies e.g. MongoDB, DynamoDB could be added. 

Notably, the support for SQL data-stores has been implemented preserving the core Nextflow capabilities to allow portable and self-contained pipeline scripts that can be developed locally, tested through CI services, and deployed at scale into production environments. 

If you have any questions or suggestions, please feel free to comment in the project discussion group at [this link](https://github.com/nextflow-io/nf-sqldb/discussions).

Credits to [Francesco Strozzi](https://twitter.com/fstrozzi) & [Raoul J.P. Bonnal](https://twitter.com/bonnalr) for having contributed to this work 🙏.
