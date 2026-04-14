import os
from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker, AsyncSession

from .base import Base
from .models import PipelineJob, StageResult

# Load from environment or default to local SQLite
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite+aiosqlite:///./data/pipeline.db")

# Ensure the data directory exists if using local SQLite
db_path = DATABASE_URL.replace("sqlite+aiosqlite:///", "")
if db_path.startswith("./"):
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

engine = create_async_engine(DATABASE_URL, echo=False)
AsyncSessionLocal = async_sessionmaker(
    bind=engine,
    class_=AsyncSession,
    expire_on_commit=False,
)

async def get_db():
    async with AsyncSessionLocal() as session:
        yield session

async def init_db():
    async with engine.begin() as conn:
        # Create all tables
        await conn.run_sync(Base.metadata.create_all)
