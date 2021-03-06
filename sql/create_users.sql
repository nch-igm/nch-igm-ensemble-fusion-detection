CREATE USER fusion_user_ro WITH PASSWORD 'PASSWORD_RO_HERE';
GRANT CONNECT ON DATABASE fusions TO fusion_user_ro;
GRANT USAGE ON SCHEMA public TO fusion_user_ro;
GRANT SELECT ON ALL TABLES IN SCHEMA public TO fusion_user_ro;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO fusion_user_ro;

CREATE USER fusion_user_rw WITH PASSWORD 'PASSWORD_RW_HERE';
GRANT CONNECT ON DATABASE fusions TO fusion_user_rw;
GRANT USAGE, SELECT ON ALL SEQUENCES IN SCHEMA public TO fusion_user_rw;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO fusion_user_rw;
